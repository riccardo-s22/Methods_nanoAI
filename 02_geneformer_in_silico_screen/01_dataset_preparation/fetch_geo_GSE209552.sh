#!/bin/bash
# Fetch the public input of 01_dataset_preparation: GEO series GSE209552
# (Garza et al. 2023), about 0.9 GB.
#
# NOT part of the original run. The original pipeline started from a manual
# download; this helper was added so the stage can start from a fresh clone.
# It downloads and unpacks, and writes the sample manifest that
# build_geneformer_h5ad_from_10x.py reads. It does no analysis.
#
#   bash 02_geneformer_in_silico_screen/01_dataset_preparation/fetch_geo_GSE209552.sh
#
# Writes (locations from ../gf_env.sh):
#   $GF_GEO/GSE209552_RAW.tar            the series supplementary archive
#   $GF_GEO/GSE209552_RAW/               its contents (10x triplets per GSM)
#   $GF_GEO/GSE209552_family.soft        the series SOFT file
#   $GF_WORK/reports/snrna_sample_manifest.txt
#       one line per snRNA-seq sample: GSM|matrix|barcodes|features
#       (the 17 samples that have a 10x matrix; TE matrices are skipped)
#
# Not fetched: the SRA Run Selector table. parse_geo_series_metadata.py merges it
# when $GF_GEO/SraRunTable.csv exists; without it the per-sample table lacks only
# the five *_sra columns. To add it, open the series' SRA study in the Run
# Selector, download the metadata table and save it under that name.
#
# Files already present are kept, so the script can be re-run.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../gf_env.sh"

GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209552"
RAW="${GF_GEO}/GSE209552_RAW"
MANIFEST="${GF_WORK}/reports/snrna_sample_manifest.txt"

fetch() {  # fetch <url> <destination>; kept if already there
    if [ -s "$2" ]; then
        echo "have    $2"
        return 0
    fi
    echo "fetch   $1"
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 -o "$2.part" "$1"
    else
        wget -O "$2.part" "$1"
    fi
    mv "$2.part" "$2"
}

mkdir -p "${RAW}" "$(dirname "${MANIFEST}")"

fetch "${GEO_URL}/suppl/GSE209552_RAW.tar" "${GF_GEO}/GSE209552_RAW.tar"
if [ ! -e "${GF_GEO}/GSE209552_family.soft" ]; then
    fetch "${GEO_URL}/soft/GSE209552_family.soft.gz" "${GF_GEO}/GSE209552_family.soft.gz"
    gunzip -k "${GF_GEO}/GSE209552_family.soft.gz"
fi

echo "unpack  ${RAW}"
tar -xf "${GF_GEO}/GSE209552_RAW.tar" -C "${RAW}" --skip-old-files

if [ -s "${MANIFEST}" ]; then
    echo "have    ${MANIFEST}"
else
    ls "${RAW}" | grep -E '_matrix\.mtx\.gz$' | grep -v 'TE_' | sort | while read -r m; do
        s="${m%_matrix.mtx.gz}"
        echo "${m%%_*}|${m}|${s}_barcodes.tsv.gz|${s}_features.tsv.gz"
    done > "${MANIFEST}"
    echo "wrote   ${MANIFEST}"
fi
echo "samples in manifest: $(wc -l < "${MANIFEST}")"
