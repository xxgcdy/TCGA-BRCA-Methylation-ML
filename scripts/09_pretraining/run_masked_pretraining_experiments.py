#!/usr/bin/env python3
"""Masked methylation pretraining, position ablation, and low-label fine-tuning.

Pretraining and label subsampling use training patients only. Validation labels
select checkpoints and thresholds; the held-out test split is never used for
model selection.
"""

import json
import math
import os
import random
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, average_precision_score, confusion_matrix, roc_auc_score, roc_curve
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
DATA = ROOT / "processed/transformer/tcga_brca_hm450_top1000.h5"
META = ROOT / "processed/transformer/sample_metadata.tsv"
OUT = ROOT / "results/masked_pretraining"
CKPT = OUT / "checkpoints"
OUT.mkdir(parents=True, exist_ok=True); CKPT.mkdir(parents=True, exist_ok=True)

DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if os.environ.get("REQUIRE_CUDA", "1") == "1" and DEVICE.type != "cuda":
    raise RuntimeError("CUDA is required for this experiment")

SEED = 20261028
PATCH_SIZE = 10
D_MODEL = 64
POSITION_MODES = os.environ.get("POSITION_MODES", "none,sequence,genomic").split(",")
LABEL_FRACTIONS = [float(x) for x in os.environ.get("LABEL_FRACTIONS", "0.10,0.25,0.50,1.00").split(",")]
FINETUNE_SEEDS = [int(x) for x in os.environ.get("FINETUNE_SEEDS", "20261027,20261028,20261029").split(",")]
PRETRAIN_EPOCHS = int(os.environ.get("PRETRAIN_EPOCHS", "50"))
FINETUNE_EPOCHS = int(os.environ.get("FINETUNE_EPOCHS", "50"))
MASK_RATE = float(os.environ.get("MASK_RATE", "0.30"))


def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    if torch.cuda.is_available(): torch.cuda.manual_seed_all(seed)


set_seed(SEED)
with h5py.File(DATA, "r") as h:
    beta = h["beta"][:].astype(np.float32)
    labels = h["label"][:].astype(np.float32)
    split = h["split_code"][:].astype(np.int8)
    mean = h["train_mean"][:].astype(np.float32)
    std = h["train_std"][:].astype(np.float32)
    chromosomes = h["chromosome"][:].astype(np.int64)
    positions = h["position"][:].astype(np.int64)
meta = pd.read_csv(META, sep="\t")

beta = np.clip((beta - mean) / std, -6, 6).astype(np.float32)


def make_patches(values, chrom, pos):
    blocks, member_masks, patch_chr, patch_mid, members_all = [], [], [], [], []
    for c in sorted(np.unique(chrom)):
        idx = np.where(chrom == c)[0]
        for start in range(0, len(idx), PATCH_SIZE):
            members = idx[start:start + PATCH_SIZE]
            block = np.zeros((len(values), PATCH_SIZE), dtype=np.float32)
            member_mask = np.zeros(PATCH_SIZE, dtype=bool)
            block[:, :len(members)] = values[:, members]
            member_mask[:len(members)] = True
            blocks.append(block); member_masks.append(member_mask); patch_chr.append(c)
            patch_mid.append(float(np.mean(pos[members]))); members_all.append(members)
    assert all(len(set(chrom[m])) == 1 for m in members_all)
    return np.stack(blocks, axis=1), np.stack(member_masks), np.asarray(patch_chr), np.asarray(patch_mid)


x, member_mask, patch_chr, patch_mid = make_patches(beta, chromosomes, positions)
N_PATCHES = x.shape[1]
member_mask_t = torch.from_numpy(member_mask).to(DEVICE)
patch_chr_t = torch.from_numpy(patch_chr).long().to(DEVICE)

# Coordinate features are fixed and independent of labels.
coord = patch_mid / 250_000_000.0
coord_features = np.stack([
    coord,
    np.sin(2 * np.pi * coord), np.cos(2 * np.pi * coord),
    np.sin(4 * np.pi * coord), np.cos(4 * np.pi * coord),
    np.sin(8 * np.pi * coord), np.cos(8 * np.pi * coord),
    np.log1p(patch_mid) / np.log(250_000_001.0),
], axis=1).astype(np.float32)
coord_features_t = torch.from_numpy(coord_features).to(DEVICE)


class MethylationEncoder(nn.Module):
    def __init__(self, position_mode):
        super().__init__()
        self.position_mode = position_mode
        self.patch = nn.Sequential(nn.Linear(PATCH_SIZE, D_MODEL), nn.LayerNorm(D_MODEL), nn.GELU())
        self.cls = nn.Parameter(torch.zeros(1, 1, D_MODEL))
        self.mask_token = nn.Parameter(torch.zeros(1, 1, D_MODEL))
        self.chrom = nn.Embedding(23, D_MODEL)
        self.sequence_pos = nn.Parameter(torch.zeros(1, N_PATCHES + 1, D_MODEL))
        self.coord_projection = nn.Sequential(nn.Linear(coord_features.shape[1], D_MODEL), nn.LayerNorm(D_MODEL))
        layer = nn.TransformerEncoderLayer(D_MODEL, 4, D_MODEL * 2, 0.15,
                                           activation="gelu", batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, 2)
        self.norm = nn.LayerNorm(D_MODEL)
        self.decoder = nn.Linear(D_MODEL, PATCH_SIZE)
        self.classifier = nn.Sequential(nn.Linear(D_MODEL, 32), nn.GELU(), nn.Dropout(0.20), nn.Linear(32, 1))
        nn.init.trunc_normal_(self.cls, std=0.02); nn.init.trunc_normal_(self.mask_token, std=0.02)
        nn.init.trunc_normal_(self.sequence_pos, std=0.02)

    def position_embedding(self, batch_size):
        zero_cls = torch.zeros(1, 1, D_MODEL, device=DEVICE)
        if self.position_mode == "none":
            return zero_cls.expand(batch_size, N_PATCHES + 1, D_MODEL)
        if self.position_mode == "sequence":
            patch_pos = self.sequence_pos[:, 1:] + self.chrom(patch_chr_t)[None]
            return torch.cat([self.sequence_pos[:, :1], patch_pos], dim=1).expand(batch_size, -1, -1)
        if self.position_mode == "genomic":
            patch_pos = self.chrom(patch_chr_t) + self.coord_projection(coord_features_t)
            return torch.cat([zero_cls, patch_pos[None]], dim=1).expand(batch_size, -1, -1)
        raise ValueError(self.position_mode)

    def encode(self, z, masked_patches=None):
        token = self.patch(z)
        if masked_patches is not None:
            token = torch.where(masked_patches[:, :, None], self.mask_token.expand_as(token), token)
        seq = torch.cat([self.cls.expand(z.shape[0], -1, -1), token], dim=1)
        return self.encoder(seq + self.position_embedding(z.shape[0]))

    def pretrain(self, z, masked_patches):
        return self.decoder(self.encode(z, masked_patches)[:, 1:])

    def classify(self, z):
        return self.classifier(self.norm(self.encode(z)[:, 0])).squeeze(1)


def loader(indices, batch_size, shuffle, seed):
    ds = TensorDataset(torch.from_numpy(x[indices]), torch.from_numpy(labels[indices]))
    gen = torch.Generator().manual_seed(seed)
    return DataLoader(ds, batch_size=batch_size, shuffle=shuffle, generator=gen, num_workers=0, pin_memory=True)


def pretrain(position_mode):
    set_seed(SEED)
    model = MethylationEncoder(position_mode).to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=5e-4, weight_decay=1e-3)
    train_idx = np.where(split == 0)[0]
    train_loader = loader(train_idx, 64, True, SEED)
    history = []
    for epoch in range(1, PRETRAIN_EPOCHS + 1):
        model.train(); losses = []
        for xb, _ in train_loader:
            xb = xb.to(DEVICE, non_blocking=True)
            masked = torch.rand((len(xb), N_PATCHES), device=DEVICE) < MASK_RATE
            masked[:, torch.randint(N_PATCHES, (1,), device=DEVICE)] = True
            prediction = model.pretrain(xb, masked)
            weight = masked[:, :, None] & member_mask_t[None]
            loss = ((prediction - xb).square() * weight).sum() / weight.sum().clamp_min(1)
            optimizer.zero_grad(set_to_none=True); loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            losses.append(float(loss.detach()))
        history.append({"position_mode": position_mode, "epoch": epoch, "masked_mse": float(np.mean(losses))})
        if epoch == 1 or epoch % 10 == 0:
            print(json.dumps(history[-1]), flush=True)
    path = CKPT / f"pretrained_{position_mode}.pt"
    torch.save({"model_state": model.state_dict(), "position_mode": position_mode,
                "mask_rate": MASK_RATE, "epochs": PRETRAIN_EPOCHS}, path)
    return path, history


def patient_subset(fraction, seed):
    train_rows = np.where(split == 0)[0]
    train_meta = meta.iloc[train_rows].copy()
    case_stratum = train_meta.groupby("case_id")["label"].agg(
        lambda z: "paired" if set(z.astype(int)) == {0, 1} else ("tumor_only" if int(z.max()) == 1 else "normal_only"))
    rng = np.random.default_rng(seed)
    selected = []
    for _, cases in case_stratum.groupby(case_stratum).groups.items():
        cases = np.asarray(list(cases), dtype=object)
        n = len(cases) if fraction >= 1 else max(1, int(round(len(cases) * fraction)))
        selected.extend(rng.choice(cases, n, replace=False).tolist())
    chosen = train_meta["case_id"].isin(selected).to_numpy()
    return train_rows[chosen]


@torch.no_grad()
def predict(model, data_loader):
    model.eval(); yy, pp = [], []
    for xb, yb in data_loader:
        prob = torch.sigmoid(model.classify(xb.to(DEVICE, non_blocking=True))).cpu().numpy()
        yy.append(yb.numpy()); pp.append(prob)
    return np.concatenate(yy).astype(int), np.concatenate(pp)


def metric_row(position_mode, initialization, fraction, seed, split_name, yy, pp, threshold, n_train):
    pred = (pp >= threshold).astype(int); tn, fp, fn, tp = confusion_matrix(yy, pred, labels=[0, 1]).ravel()
    return {"position_mode": position_mode, "initialization": initialization, "label_fraction": fraction,
            "seed": seed, "split": split_name, "n_train": n_train, "n": len(yy), "threshold": threshold,
            "roc_auc": roc_auc_score(yy, pp), "average_precision": average_precision_score(yy, pp),
            "accuracy": accuracy_score(yy, pred), "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp)}


def finetune(position_mode, initialization, fraction, seed, pretrained_path=None):
    set_seed(seed)
    model = MethylationEncoder(position_mode).to(DEVICE)
    if pretrained_path is not None:
        state = torch.load(pretrained_path, map_location=DEVICE, weights_only=True)["model_state"]
        model.load_state_dict(state)
        model.classifier.apply(lambda module: module.reset_parameters() if hasattr(module, "reset_parameters") else None)
    train_idx = patient_subset(fraction, seed)
    val_idx = np.where(split == 1)[0]; test_idx = np.where(split == 2)[0]
    train_loader = loader(train_idx, 64, True, seed)
    val_loader = loader(val_idx, 128, False, seed); test_loader = loader(test_idx, 128, False, seed)
    y_train = labels[train_idx].astype(int)
    pos_weight = torch.tensor([(y_train == 0).sum() / max((y_train == 1).sum(), 1)], device=DEVICE)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-3)
    best_auc, stale, best_epoch = -np.inf, 0, 0
    path = CKPT / f"ft_{position_mode}_{initialization}_f{fraction:.2f}_s{seed}.pt"
    history = []
    for epoch in range(1, FINETUNE_EPOCHS + 1):
        model.train(); losses = []
        for xb, yb in train_loader:
            xb = xb.to(DEVICE, non_blocking=True); yb = yb.to(DEVICE, non_blocking=True)
            loss = criterion(model.classify(xb), yb)
            optimizer.zero_grad(set_to_none=True); loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step(); losses.append(float(loss.detach()))
        vy, vp = predict(model, val_loader); val_auc = roc_auc_score(vy, vp)
        history.append({"position_mode": position_mode, "initialization": initialization,
                        "label_fraction": fraction, "seed": seed, "epoch": epoch,
                        "train_loss": float(np.mean(losses)), "val_auc": val_auc})
        if val_auc > best_auc + 1e-5:
            best_auc, stale, best_epoch = val_auc, 0, epoch
            torch.save(model.state_dict(), path)
        else:
            stale += 1
            if stale >= 8: break
    model.load_state_dict(torch.load(path, map_location=DEVICE, weights_only=True))
    vy, vp = predict(model, val_loader); ty, tp = predict(model, test_loader)
    fpr, tpr, thresholds = roc_curve(vy, vp); threshold = float(thresholds[np.argmax(tpr - fpr)])
    rows = [metric_row(position_mode, initialization, fraction, seed, "validation", vy, vp, threshold, len(train_idx)),
            metric_row(position_mode, initialization, fraction, seed, "test", ty, tp, threshold, len(train_idx))]
    predictions = pd.concat([
        meta.iloc[val_idx][["sample_id", "case_id", "label"]].assign(
            position_mode=position_mode, initialization=initialization, label_fraction=fraction, seed=seed,
            split="validation", probability=vp),
        meta.iloc[test_idx][["sample_id", "case_id", "label"]].assign(
            position_mode=position_mode, initialization=initialization, label_fraction=fraction, seed=seed,
            split="test", probability=tp),
    ])
    return rows, history, predictions, best_epoch


started = time.time(); pretrain_history = []; checkpoints = {}
for mode in POSITION_MODES:
    checkpoints[mode], hist = pretrain(mode); pretrain_history.extend(hist)

metric_rows, finetune_history, prediction_rows, best_epochs = [], [], [], []
for mode in POSITION_MODES:
    for initialization in ["random", "pretrained"]:
        for fraction in LABEL_FRACTIONS:
            for seed in FINETUNE_SEEDS:
                rows, hist, pred, best_epoch = finetune(
                    mode, initialization, fraction, seed,
                    checkpoints[mode] if initialization == "pretrained" else None)
                metric_rows.extend(rows); finetune_history.extend(hist); prediction_rows.append(pred)
                best_epochs.append({"position_mode": mode, "initialization": initialization,
                                    "label_fraction": fraction, "seed": seed, "best_epoch": best_epoch})
                print(json.dumps(rows[-1]), flush=True)

metrics_df = pd.DataFrame(metric_rows)
metrics_df.to_csv(OUT / "finetune_metrics_all_runs.tsv", sep="\t", index=False)
pd.DataFrame(pretrain_history).to_csv(OUT / "pretraining_history.tsv", sep="\t", index=False)
pd.DataFrame(finetune_history).to_csv(OUT / "finetune_history.tsv.gz", sep="\t", index=False, compression="gzip")
pd.concat(prediction_rows, ignore_index=True).to_csv(OUT / "predictions_all_runs.tsv.gz", sep="\t", index=False, compression="gzip")
pd.DataFrame(best_epochs).to_csv(OUT / "best_epochs.tsv", sep="\t", index=False)

test = metrics_df.query("split == 'test'")
group_cols = ["position_mode", "initialization", "label_fraction"]
summary_table = test.groupby(group_cols).agg(
    n_seeds=("seed", "nunique"), roc_auc_mean=("roc_auc", "mean"), roc_auc_sd=("roc_auc", "std"),
    accuracy_mean=("accuracy", "mean"), sensitivity_mean=("sensitivity", "mean"),
    specificity_mean=("specificity", "mean")).reset_index()
summary_table.to_csv(OUT / "test_summary_mean_sd.tsv", sep="\t", index=False)

run_summary = {
    "device": str(DEVICE), "torch_version": torch.__version__, "cuda_version": torch.version.cuda,
    "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
    "gpu_capability": torch.cuda.get_device_capability(0) if torch.cuda.is_available() else None,
    "n_samples": int(len(x)), "n_patches": int(N_PATCHES), "cross_chromosome_patches": 0,
    "position_modes": POSITION_MODES, "label_fractions": LABEL_FRACTIONS,
    "finetune_seeds": FINETUNE_SEEDS, "pretrain_epochs": PRETRAIN_EPOCHS,
    "finetune_max_epochs": FINETUNE_EPOCHS, "mask_rate": MASK_RATE,
    "elapsed_minutes": (time.time() - started) / 60,
}
(OUT / "run_summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
print(summary_table.to_string(index=False), flush=True)
print(json.dumps(run_summary, indent=2), flush=True)
