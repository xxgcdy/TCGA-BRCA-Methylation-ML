#!/usr/bin/env python3
import json
import math
import random
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import (accuracy_score, average_precision_score,
                             confusion_matrix, roc_auc_score, roc_curve)
from torch import nn
from torch.utils.data import DataLoader, TensorDataset, WeightedRandomSampler

SEED = 20261027
random.seed(SEED); np.random.seed(SEED); torch.manual_seed(SEED)
torch.set_num_threads(8); torch.set_num_interop_threads(1)
ROOT = Path(os.environ.get("METHYLATION_ROOT", Path.cwd()))
DATA = ROOT / "processed/transformer/tcga_brca_hm450_top1000.h5"
OUT = ROOT / "results/transformer"
OUT.mkdir(parents=True, exist_ok=True)

with h5py.File(DATA, "r") as h:
    x = h["beta"][:].astype(np.float32)
    y = h["label"][:].astype(np.float32)
    split = h["split_code"][:]
    mean = h["train_mean"][:].astype(np.float32)
    std = h["train_std"][:].astype(np.float32)
    chrom = h["chromosome"][:].astype(np.int64)

x = np.clip((x - mean) / std, -6, 6).astype(np.float32)
patch_size = 10
n_patches = x.shape[1] // patch_size
x = x[:, :n_patches * patch_size].reshape(len(x), n_patches, patch_size)
patch_chrom = chrom[:n_patches * patch_size].reshape(n_patches, patch_size)[:, 0]

class PatchTransformer(nn.Module):
    def __init__(self, n_patches, patch_size, d_model=64, heads=4, layers=2, dropout=0.20):
        super().__init__()
        self.patch = nn.Sequential(nn.Linear(patch_size, d_model), nn.LayerNorm(d_model), nn.GELU())
        self.cls = nn.Parameter(torch.zeros(1, 1, d_model))
        self.pos = nn.Parameter(torch.zeros(1, n_patches + 1, d_model))
        self.chrom = nn.Embedding(23, d_model)
        enc_layer = nn.TransformerEncoderLayer(d_model, heads, d_model * 2, dropout,
                                               activation="gelu", batch_first=True, norm_first=True)
        self.encoder = nn.TransformerEncoder(enc_layer, layers)
        self.norm = nn.LayerNorm(d_model)
        self.head = nn.Sequential(nn.Linear(d_model, 32), nn.GELU(), nn.Dropout(dropout), nn.Linear(32, 1))
        nn.init.trunc_normal_(self.cls, std=0.02); nn.init.trunc_normal_(self.pos, std=0.02)

    def forward(self, z, patch_chr):
        b = z.shape[0]
        tok = self.patch(z) + self.chrom(patch_chr)[None, :, :]
        cls = self.cls.expand(b, -1, -1)
        z = torch.cat([cls, tok], dim=1) + self.pos
        z = self.encoder(z)
        return self.head(self.norm(z[:, 0])).squeeze(1)

def loader_for(code, train=False):
    idx = np.where(split == code)[0]
    ds = TensorDataset(torch.from_numpy(x[idx]), torch.from_numpy(y[idx]))
    if train:
        counts = np.bincount(y[idx].astype(int), minlength=2)
        weights = np.where(y[idx] == 1, 1 / counts[1], 1 / counts[0])
        sampler = WeightedRandomSampler(torch.as_tensor(weights, dtype=torch.double), len(idx), replacement=True)
        return DataLoader(ds, batch_size=32, sampler=sampler, num_workers=0), idx
    return DataLoader(ds, batch_size=64, shuffle=False, num_workers=0), idx

train_loader, train_idx = loader_for(0, True)
val_loader, val_idx = loader_for(1)
test_loader, test_idx = loader_for(2)
patch_chr_t = torch.from_numpy(patch_chrom).long()
model = PatchTransformer(n_patches, patch_size)
optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4, weight_decay=1e-3)
criterion = nn.BCEWithLogitsLoss()
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="max", factor=0.5, patience=3)

@torch.no_grad()
def predict(loader):
    model.eval(); probs = []; labels = []
    for xb, yb in loader:
        probs.append(torch.sigmoid(model(xb, patch_chr_t)).cpu().numpy()); labels.append(yb.numpy())
    return np.concatenate(labels).astype(int), np.concatenate(probs)

history = []
best_auc = -np.inf; best_epoch = 0; patience = 10; stale = 0
for epoch in range(1, 51):
    model.train(); losses = []
    for xb, yb in train_loader:
        optimizer.zero_grad(set_to_none=True)
        loss = criterion(model(xb, patch_chr_t), yb)
        loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 1.0); optimizer.step()
        losses.append(float(loss.detach()))
    vy, vp = predict(val_loader)
    va = roc_auc_score(vy, vp)
    scheduler.step(va)
    history.append({"epoch": epoch, "train_loss": float(np.mean(losses)), "val_auc": float(va),
                    "lr": optimizer.param_groups[0]["lr"]})
    print(json.dumps(history[-1]), flush=True)
    if va > best_auc + 1e-5:
        best_auc = va; best_epoch = epoch; stale = 0
        torch.save({"model_state": model.state_dict(), "architecture": {"patch_size": patch_size,
                   "n_patches": n_patches, "d_model": 64, "heads": 4, "layers": 2},
                   "seed": SEED, "best_epoch": best_epoch, "best_val_auc": best_auc},
                   OUT / "best_model.pt")
    else:
        stale += 1
        if stale >= patience: break

checkpoint = torch.load(OUT / "best_model.pt", map_location="cpu", weights_only=False)
model.load_state_dict(checkpoint["model_state"])
vy, vp = predict(val_loader); ty, tp = predict(test_loader)
fpr, tpr, thresholds = roc_curve(vy, vp)
threshold = float(thresholds[np.argmax(tpr - fpr)])

def metric_row(name, yy, pp):
    pred = (pp >= threshold).astype(int)
    tn, fp, fn, tpv = confusion_matrix(yy, pred, labels=[0, 1]).ravel()
    return {"split": name, "n": len(yy), "n_tumor": int(yy.sum()), "n_normal": int((yy == 0).sum()),
            "threshold": threshold, "roc_auc": float(roc_auc_score(yy, pp)),
            "average_precision": float(average_precision_score(yy, pp)),
            "accuracy": float(accuracy_score(yy, pred)), "sensitivity": float(tpv / (tpv + fn)),
            "specificity": float(tn / (tn + fp))}

metrics = [metric_row("validation", vy, vp), metric_row("test", ty, tp)]
pd.DataFrame(metrics).to_csv(OUT / "transformer_metrics.tsv", sep="\t", index=False)
sample_meta = pd.read_csv(ROOT / "processed/transformer/sample_metadata.tsv", sep="\t")
preds = pd.concat([
    sample_meta.iloc[val_idx].assign(model="CpGPatchTransformer", probability=vp),
    sample_meta.iloc[test_idx].assign(model="CpGPatchTransformer", probability=tp),
])
preds.to_csv(OUT / "transformer_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
pd.DataFrame(history).to_csv(OUT / "training_history.tsv", sep="\t", index=False)
summary = {"seed": SEED, "best_epoch": best_epoch, "best_val_auc": best_auc,
           "architecture": checkpoint["architecture"], "n_parameters": sum(p.numel() for p in model.parameters()),
           "metrics": metrics}
(OUT / "transformer_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)

