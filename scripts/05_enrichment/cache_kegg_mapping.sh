#!/usr/bin/env bash
set -euo pipefail

root="${METHYLATION_ROOT:-$(pwd)}"
out="$root/annotation/kegg"
mkdir -p "$out"

curl -L --fail --retry 5 --retry-delay 5 --connect-timeout 30 \
  --speed-time 120 --speed-limit 100 \
  -o "$out/hsa_pathway_links.tsv" \
  https://rest.kegg.jp/link/pathway/hsa

curl -L --fail --retry 5 --retry-delay 5 --connect-timeout 30 \
  --speed-time 120 --speed-limit 100 \
  -o "$out/hsa_pathways.tsv" \
  https://rest.kegg.jp/list/pathway/hsa

wc -l "$out/hsa_pathway_links.tsv" "$out/hsa_pathways.tsv"
