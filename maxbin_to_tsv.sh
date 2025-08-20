#!/bin/bash
for bin in "$1"; do
  bin_name=$(basename "$bin" .fasta)
  awk -v bname="$bin_name" '/^>/ {gsub(/^>/,""); print $1"\t"bname}' "$bin"
done