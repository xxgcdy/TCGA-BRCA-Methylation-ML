#!/usr/bin/env python3
"""Train chromosome-safe CpG-patch Transformers across multiple seeds."""

import json
import os
import random
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, average_precision_score, confusion_matrix, roc_auc_score, roc_curve
from torch import nn
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
DATA = ROOT / "processed/transformer/tcga_brca_hm450_top1000.h5"
OUT = ROOT / "results/transformer_chromosome_safe"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = [int(x) for x in os.environ.get("TRANSFORMER_SEEDS", "20261027,20261028,20261029,20261030,20261031").split(",")]
PATCH_SIZE = 10
torch.set_num_threads(8)
torch.set_num_interop_threads(1)

with h5py.File(DATA, "r") as h:
    raw_x = h["beta"][:].astype(np.float32)
    y = h["label"][:].astype(np.float32)
    split = h["split_code"][:]
    mean = h["train_mean"][:].astype(np.float32)
    std = h["train_std"][:].astype(np.float32)
    chrom = h["chromosome"][:].astype(np.int64)
    position = h["position"][:].astype(np.int64)

raw_x = np.clip((raw_x - mean) / std, -6, 6).astype(np.float32)


def chromosome_safe_patches(values, chromosomes, positions, patch_size=10):
    chunks, masks, patch_chr, start_pos, end_pos, member_indices = [], [], [], [], [], []
    for chromosome in sorted(np.unique(chromosomes)):
        idx = np.where(chromosomes == chromosome)[0]
        for start in range(0, len(idx), patch_size):
            members = idx[start:start + patch_size]
            block = np.zeros((len(values), patch_size), dtype=np.float32)
            mask = np.zeros(patch_size, dtype=bool)
            block[:, :len(members)] = values[:, members]
            mask[:len(members)] = True
            chunks.append(block); masks.append(mask); patch_chr.append(chromosome)
            start_pos.append(int(positions[members[0]])); end_pos.append(int(positions[members[-1]]))
            member_indices.append(members.tolist())
    return (np.stack(chunks, axis=1), np.stack(masks), np.asarray(patch_chr, dtype=np.int64),
            np.asarray(start_pos), np.asarray(end_pos), member_indices)


x, member_mask, patch_chr, patch_start, patch_end, member_indices = chromosome_safe_patches(
    raw_x, chrom, position, PATCH_SIZE)
assert all(len(set(chrom[m])) == 1 for m in member_indices)

patch_table = pd.DataFrame({"patch_index": np.arange(len(patch_chr)), "chromosome": patch_chr,
                            "start_hg19": patch_start, "end_hg19": patch_end,
                            "n_cpg": member_mask.sum(axis=1)})
patch_table.to_csv(OUT / "patch_metadata.tsv", sep="\t", index=False)


class PatchTransformer(nn.Module):
    def __init__(self, n_patches, patch_size, d_model=64, heads=4, layers=2, dropout=0.20):
        super().__init__()
        self.patch = nn.Sequential(nn.Linear(patch_size, d_model), nn.LayerNorm(d_model), nn.GELU())
        self.cls = nn.Parameter(torch.zeros(1, 1, d_model))
        self.pos = nn.Parameter(torch.zeros(1, n_patches + 1, d_model))
        self.chrom = nn.Embedding(23, d_model)
        layer = nn.TransformerEncoderLayer(d_model, heads, d_model * 2, dropout,
                                           activation="gelu", batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, layers)
        self.norm = nn.LayerNorm(d_model)
        self.head = nn.Sequential(nn.Linear(d_model, 32), nn.GELU(), nn.Dropout(dropout), nn.Linear(32, 1))
        nn.init.trunc_normal_(self.cls, std=0.02); nn.init.trunc_normal_(self.pos, std=0.02)

    def forward(self, z, chromosomes):
        token = self.patch(z) + self.chrom(chromosomes)[None, :, :]
        z = torch.cat([self.cls.expand(z.shape[0], -1, -1), token], dim=1) + self.pos
        return self.head(self.norm(self.encoder(z)[:, 0])).squeeze(1)


def make_loader(code, seed, train=False):
    idx = np.where(split == code)[0]
    ds = TensorDataset(torch.from_numpy(x[idx]), torch.from_numpy(y[idx]))
    if train:
        counts = np.bincount(y[idx].astype(int), minlength=2)
        weights = np.where(y[idx] == 1, 1 / counts[1], 1 / counts[0])
        gen = torch.Generator().manual_seed(seed)
        sampler = WeightedRandomSampler(torch.as_tensor(weights, dtype=torch.double), len(idx), True, generator=gen)
        return DataLoader(ds, batch_size=32, sampler=sampler, num_workers=0), idx
    return DataLoader(ds, batch_size=64, shuffle=False, num_workers=0), idx


def one_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    train_loader, _ = make_loader(0, seed, True)
    val_loader, val_idx = make_loader(1, seed)
    test_loader, test_idx = make_loader(2, seed)
    patch_chr_t = torch.from_numpy(patch_chr).long()
    model = PatchTransformer(x.shape[1], PATCH_SIZE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-3)
    criterion = nn.BCEWithLogitsLoss()
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="max", factor=0.5, patience=3)

    @torch.no_grad()
    def predict(loader):
        model.eval(); probs, labels = [], []
        for xb, yb in loader:
            probs.append(torch.sigmoid(model(xb, patch_chr_t)).numpy()); labels.append(yb.numpy())
        return np.concatenate(labels).astype(int), np.concatenate(probs)

    history, best_auc, best_epoch, stale = [], -np.inf, 0, 0
    checkpoint = OUT / f"seed_{seed}_best.pt"
    for epoch in range(1, 51):
        model.train(); losses = []
        for xb, yb in train_loader:
            optimizer.zero_grad(set_to_none=True)
            loss = criterion(model(xb, patch_chr_t), yb)
            loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
            losses.append(float(loss.detach()))
        vy, vp = predict(val_loader); val_auc = roc_auc_score(vy, vp); scheduler.step(val_auc)
        history.append({"seed": seed, "epoch": epoch, "train_loss": np.mean(losses), "val_auc": val_auc,
                        "lr": optimizer.param_groups[0]["lr"]})
        if val_auc > best_auc + 1e-5:
            best_auc, best_epoch, stale = val_auc, epoch, 0
            torch.save(model.state_dict(), checkpoint)
        else:
            stale += 1
            if stale >= 10: break
    model.load_state_dict(torch.load(checkpoint, map_location="cpu", weights_only=True))
    vy, vp = predict(val_loader); ty, tp = predict(test_loader)
    threshold = float(roc_curve(vy, vp)[2][np.argmax(roc_curve(vy, vp)[1] - roc_curve(vy, vp)[0])])
    return model, history, val_idx, test_idx, vy, vp, ty, tp, threshold, best_epoch


def metrics(seed, split_name, yy, pp, threshold):
    pred = (pp >= threshold).astype(int); tn, fp, fn, tp = confusion_matrix(yy, pred, labels=[0, 1]).ravel()
    return {"seed": seed, "split": split_name, "n": len(yy), "threshold": threshold,
            "roc_auc": roc_auc_score(yy, pp), "average_precision": average_precision_score(yy, pp),
            "accuracy": accuracy_score(yy, pred), "sensitivity": tp / (tp + fn), "specificity": tn / (tn + fp)}


all_metrics, all_history, val_prob, test_prob, best_epochs = [], [], [], [], []
for seed in SEEDS:
    model, history, val_idx, test_idx, vy, vp, ty, tp, threshold, best_epoch = one_seed(seed)
    all_history.extend(history); val_prob.append(vp); test_prob.append(tp); best_epochs.append(best_epoch)
    all_metrics += [metrics(seed, "validation", vy, vp, threshold), metrics(seed, "test", ty, tp, threshold)]

val_ensemble = np.mean(val_prob, axis=0); test_ensemble = np.mean(test_prob, axis=0)
fpr, tpr, thresholds = roc_curve(vy, val_ensemble)
ensemble_threshold = float(thresholds[np.argmax(tpr - fpr)])
all_metrics += [metrics("ensemble", "validation", vy, val_ensemble, ensemble_threshold),
                metrics("ensemble", "test", ty, test_ensemble, ensemble_threshold)]

pd.DataFrame(all_metrics).to_csv(OUT / "metrics_all_seeds.tsv", sep="\t", index=False)
pd.DataFrame(all_history).to_csv(OUT / "training_history_all_seeds.tsv", sep="\t", index=False)
sample_meta = pd.read_csv(ROOT / "processed/transformer/sample_metadata.tsv", sep="\t")
pred = pd.concat([sample_meta.iloc[val_idx].assign(model="ChromosomeSafeTransformerEnsemble", probability=val_ensemble),
                  sample_meta.iloc[test_idx].assign(model="ChromosomeSafeTransformerEnsemble", probability=test_ensemble)])
pred.to_csv(OUT / "ensemble_predictions.tsv.gz", sep="\t", index=False, compression="gzip")

test_rows = pd.DataFrame(all_metrics).query("split == 'test' and seed != 'ensemble'")
summary = {"seeds": SEEDS, "patch_size": PATCH_SIZE, "n_patches": int(x.shape[1]),
           "cross_chromosome_patches": 0, "partial_patches": int((member_mask.sum(axis=1) < PATCH_SIZE).sum()),
           "best_epochs": best_epochs, "test_mean_across_seeds": test_rows.select_dtypes("number").mean().to_dict(),
           "test_sd_across_seeds": test_rows.select_dtypes("number").std(ddof=1).to_dict(),
           "ensemble_metrics": [x for x in all_metrics if x["seed"] == "ensemble"]}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)
