# HOMOLOGOUS PHYLOGENETIC TREES
```bash
#!/bin/bash
set -euo pipefail

# Configuration
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# Directories
ID_FILE="~/binning/filtered/all/list_no_fa.txt"


HMM_OUT_MAGS="/home/omar/Proyecto/analisis/hmmer/mags"
PHMMER_OUT_MAGS="/home/omar/Proyecto/analisis/phmmer/mags"
HOMOLOGY_HIT="/home/omar/Proyecto/analisis/hom_hit"
DOMINIOS="/home/omar/Proyecto/analisis/dominios"

# Create output folders
mkdir -p  "$HOMOLOGY_HIT" "$DOMINIOS"

# MOVE ALL _hit.faa
find "$PHMMER_OUT_MAGS/" -type f -name '*_hits.faa' -exec cp {} "$HOMOLOGY_HIT" \;
find "$HMM_OUT_MAGS/"   -type f -name '*_hits.faa' -exec cp {} "$HOMOLOGY_HIT" \;


# ========================
# ALIGN
# ========================
cd "$HOMOLOGY_HIT"

# Extract family names (second part between "_" and before "_aligned_hits.faa")
FAMILIAS_OUT="~/analisis/familias.txt"
> "$FAMILIAS_OUT"

for id in $(cat "$ID_FILE"); do
  for file in "$HOMOLOGY_HIT"/*"${id}"*.faa; do
    base=$(basename "$file" .faa)
    familia=${base#${id}_}  # quitar prefijo del bin
    echo "$familia" >> "$FAMILIAS_OUT"
  done
done

mkdir -p "$HOMOLOGY_HIT/por_familia"

# concatenate
while read -r fam; do
  grep_pattern="*_""${fam}"_aligned_hits.faa
  cat $grep_pattern > "$HOMOLOGY_HIT/por_familia/${fam}.faa"
done < "$FAMILIAS_OUT"

mkdir -p "$HOMOLOGY_HIT/alineamientos"

# align for mechanism
for f in "$HOMOLOGY_HIT/por_familia/"*.faa; do
  base=$(basename "$f" .faa)
  conda run -n align_env mafft --auto "$f" \
  > "$HOMOLOGY_HIT/alineamientos/${base}_aligned.faa"
done

# ================================
# HOMOLOGOUS PHYLOGENETIC TREES
# ================================

# clean with trimAL
mkdir -p "$HOMOLOGY_HIT/alineamientos_trimmed"

for aln in "$HOMOLOGY_HIT/alineamientos/"*.faa; do
  base=$(basename "$aln" _aligned.faa)

  conda run -n Hom_validation trimal \
    -in "$aln" \
    -out "$HOMOLOGY_HIT/alineamientos_trimmed/${base}_trimmed.faa" \
    -automated1
done

# tree
mkdir -p "$HOMOLOGY_HIT/trees_hom"

for aln in "$HOMOLOGY_HIT/alineamientos_trimmed/"*.faa; do
  base=$(basename "$aln" .faa)

  conda run -n tree_env iqtree \
    -s "$aln" \
    -m LG+G+F \
    -bb 1000 \
    -alrt 1000 \
    -nt $THREADS_PER_JOB \
    -pre "$HOMOLOGY_HIT/trees_hom/${base}"
done
```
