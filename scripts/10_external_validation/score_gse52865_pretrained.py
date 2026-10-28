#!/usr/bin/env python3
"""Frozen external evaluation of full-label position/pretraining ablations."""

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score, average_precision_score, confusion_matrix, roc_auc_score, roc_curve
from torch import nn

ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
TCGA = ROOT / "processed/transformer/tcga_brca_hm450_top1000.h5"
EXT = ROOT / "external/GSE52865/processed/gse52865_top1000.h5"
TRAIN_OUT = ROOT / "results/masked_pretraining"
OUT = ROOT / "results/external_validation/GSE52865"
OUT.mkdir(parents=True, exist_ok=True)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
PATCH_SIZE, D_MODEL = 10, 64

with h5py.File(TCGA, "r") as h:
    train_mean = h["train_mean"][:].astype(np.float32)
    train_std = h["train_std"][:].astype(np.float32)
    train_median = h["train_median"][:].astype(np.float32)
    chrom = h["chromosome"][:].astype(np.int64)
    pos = h["position"][:].astype(np.int64)
with h5py.File(EXT, "r") as h:
    ext_x = h["beta"][:].astype(np.float32)
    ext_y = h["label"][:].astype(int)
    sample_ids = np.asarray([x.decode() for x in h["sample_id"][:]])

missing = np.where(np.isnan(ext_x))
ext_x[missing] = train_median[missing[1]]
ext_x = np.clip((ext_x - train_mean) / train_std, -6, 6).astype(np.float32)


def make_patches(values):
    blocks, patch_chr, patch_mid = [], [], []
    for c in sorted(np.unique(chrom)):
        idx = np.where(chrom == c)[0]
        for start in range(0, len(idx), PATCH_SIZE):
            members = idx[start:start + PATCH_SIZE]
            block = np.zeros((len(values), PATCH_SIZE), dtype=np.float32)
            block[:, :len(members)] = values[:, members]
            blocks.append(block); patch_chr.append(c); patch_mid.append(float(np.mean(pos[members])))
    return np.stack(blocks, axis=1), np.asarray(patch_chr), np.asarray(patch_mid)


x, patch_chr, patch_mid = make_patches(ext_x)
N_PATCHES = x.shape[1]
patch_chr_t = torch.from_numpy(patch_chr).long().to(DEVICE)
coord = patch_mid / 250_000_000.0
coord_features = np.stack([coord, np.sin(2*np.pi*coord), np.cos(2*np.pi*coord),
                           np.sin(4*np.pi*coord), np.cos(4*np.pi*coord),
                           np.sin(8*np.pi*coord), np.cos(8*np.pi*coord),
                           np.log1p(patch_mid) / np.log(250_000_001.0)], axis=1).astype(np.float32)
coord_features_t = torch.from_numpy(coord_features).to(DEVICE)


class MethylationEncoder(nn.Module):
    def __init__(self, position_mode):
        super().__init__(); self.position_mode = position_mode
        self.patch = nn.Sequential(nn.Linear(PATCH_SIZE, D_MODEL), nn.LayerNorm(D_MODEL), nn.GELU())
        self.cls = nn.Parameter(torch.zeros(1, 1, D_MODEL)); self.mask_token = nn.Parameter(torch.zeros(1, 1, D_MODEL))
        self.chrom = nn.Embedding(23, D_MODEL); self.sequence_pos = nn.Parameter(torch.zeros(1, N_PATCHES + 1, D_MODEL))
        self.coord_projection = nn.Sequential(nn.Linear(coord_features.shape[1], D_MODEL), nn.LayerNorm(D_MODEL))
        layer = nn.TransformerEncoderLayer(D_MODEL, 4, D_MODEL*2, 0.15, activation="gelu", batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(layer, 2); self.norm = nn.LayerNorm(D_MODEL)
        self.decoder = nn.Linear(D_MODEL, PATCH_SIZE)
        self.classifier = nn.Sequential(nn.Linear(D_MODEL, 32), nn.GELU(), nn.Dropout(0.20), nn.Linear(32, 1))

    def position_embedding(self, batch_size):
        zero = torch.zeros(1, 1, D_MODEL, device=DEVICE)
        if self.position_mode == "none": return zero.expand(batch_size, N_PATCHES + 1, D_MODEL)
        if self.position_mode == "sequence":
            z = self.sequence_pos[:, 1:] + self.chrom(patch_chr_t)[None]
            return torch.cat([self.sequence_pos[:, :1], z], 1).expand(batch_size, -1, -1)
        z = self.chrom(patch_chr_t) + self.coord_projection(coord_features_t)
        return torch.cat([zero, z[None]], 1).expand(batch_size, -1, -1)

    def classify(self, z):
        token = self.patch(z); seq = torch.cat([self.cls.expand(len(z), -1, -1), token], 1)
        encoded = self.encoder(seq + self.position_embedding(len(z)))
        return self.classifier(self.norm(encoded[:, 0])).squeeze(1)


@torch.no_grad()
def model_probability(mode, initialization, seed):
    model = MethylationEncoder(mode).to(DEVICE)
    path = TRAIN_OUT / "checkpoints" / f"ft_{mode}_{initialization}_f1.00_s{seed}.pt"
    model.load_state_dict(torch.load(path, map_location=DEVICE, weights_only=True)); model.eval()
    probs = []
    for start in range(0, len(x), 128):
        xb = torch.from_numpy(x[start:start+128]).to(DEVICE)
        probs.append(torch.sigmoid(model.classify(xb)).cpu().numpy())
    return np.concatenate(probs)


internal = pd.read_csv(TRAIN_OUT / "predictions_all_runs.tsv.gz", sep="\t")
internal = internal.query("label_fraction == 1.0 and split == 'validation'")
rows, prediction_rows = [], []
for mode in ["none", "sequence", "genomic"]:
    for initialization in ["random", "pretrained"]:
        external_probs = []
        validation = internal.query("position_mode == @mode and initialization == @initialization")
        pivot = validation.pivot(index="sample_id", columns="seed", values="probability")
        validation_label = validation.drop_duplicates("sample_id").set_index("sample_id").loc[pivot.index, "label"].to_numpy()
        validation_ensemble = pivot.mean(axis=1).to_numpy()
        fpr, tpr, thresholds = roc_curve(validation_label, validation_ensemble)
        threshold = float(thresholds[np.argmax(tpr - fpr)])
        for seed in [20261027, 20261028, 20261029]:
            external_probs.append(model_probability(mode, initialization, seed))
        probability = np.mean(external_probs, axis=0)
        pred = (probability >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(ext_y, pred, labels=[0, 1]).ravel()
        rows.append({"position_mode": mode, "initialization": initialization, "label_fraction": 1.0,
                     "ensemble_seeds": 3, "n": len(ext_y), "threshold": threshold,
                     "roc_auc": roc_auc_score(ext_y, probability),
                     "average_precision": average_precision_score(ext_y, probability),
                     "accuracy": accuracy_score(ext_y, pred), "sensitivity": tp/(tp+fn), "specificity": tn/(tn+fp)})
        prediction_rows.append(pd.DataFrame({"sample_id": sample_ids, "label": ext_y, "probability": probability,
                                             "position_mode": mode, "initialization": initialization,
                                             "threshold": threshold}))

pd.DataFrame(rows).to_csv(OUT / "pretraining_ablation_external_metrics.tsv", sep="\t", index=False)
pd.concat(prediction_rows).to_csv(OUT / "pretraining_ablation_external_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
summary = {"cohort": "GSE52865", "frozen_tcga_features": 1000, "models": rows}
(OUT / "pretraining_ablation_external_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(pd.DataFrame(rows).to_string(index=False))
