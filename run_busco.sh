#!/bin/bash
ID=$1
THREADS=$2
DASTOOL_OUT=$3
BUSCO_OUT=$4
LINEAGE=$5

BIN_DIR="${DASTOOL_OUT}/${ID}_bins/${ID}_dastool_DASTool_bins"

echo "==> BUSCO for $ID"
mkdir -p "${BUSCO_OUT}/${ID}"

for bin in "${BIN_DIR}"/*.fa; do
  bin_name=$(basename "$bin" .fa)
  busco \
    -i "$bin" \
    -o "$bin_name" \
    -l "$LINEAGE" \
    -m genome \
    --cpu "$THREADS" \
    --out_path "${BUSCO_OUT}/${ID}"
done