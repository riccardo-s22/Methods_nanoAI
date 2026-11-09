#!/bin/bash
# Organize per-GSM 10x triplets (barcodes/features/matrix) into standard folders.
set -euo pipefail

BASE=${BASE:-/mnt/d/TBI_SC/Geneformer_TBI}
RAW=${RAW:-/mnt/d/TBI_SC/GSE209552_RAW}
OUT=$BASE/10x_by_sample
REPORT=$BASE/reports/organized_10x_files.txt

mkdir -p "$OUT" "$BASE/reports"
: > "$REPORT"

# GSMs that actually have a matrix.mtx.gz (snRNA-seq samples only).
for gsm in $(ls "$RAW" | grep -oE 'GSM[0-9]+' | sort -u); do
    matrix_files=$(ls "$RAW" 2>/dev/null | grep -E "^${gsm}_.*_matrix\.mtx\.gz$" | grep -v 'TE_raw_matrix' || true)
    if [ -z "$matrix_files" ]; then
        continue
    fi
    mkdir -p "$OUT/$gsm"

    for f in "$RAW"/${gsm}_*; do
        [ -e "$f" ] || continue
        ln -sf "$f" "$OUT/$gsm/$(basename "$f")"
    done

    cd "$OUT/$gsm"
    b=$(ls ${gsm}_*barcodes.tsv.gz 2>/dev/null | grep -v 'TE_' | head -n 1 || true)
    feat=$(ls ${gsm}_*features.tsv.gz 2>/dev/null | grep -v 'TE_' | head -n 1 || true)
    g=$(ls ${gsm}_*genes.tsv.gz 2>/dev/null | grep -v 'TE_' | head -n 1 || true)
    m=$(ls ${gsm}_*matrix.mtx.gz 2>/dev/null | grep -v 'TE_' | head -n 1 || true)

    [ -n "$b" ] && ln -sf "$b" barcodes.tsv.gz
    [ -n "$feat" ] && ln -sf "$feat" features.tsv.gz
    [ -z "$feat" ] && [ -n "$g" ] && ln -sf "$g" features.tsv.gz
    [ -n "$m" ] && ln -sf "$m" matrix.mtx.gz
    cd "$BASE"
done

find "$OUT" -maxdepth 2 \( -type l -o -type f \) | sort > "$REPORT"
echo "Wrote $REPORT"
