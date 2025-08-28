# HOMOLOGOUS PHYLOGENETIC TREES
```bash
#!/bin/bash
set -euo pipefail


# ---------------------------
# CONFIGURATION
# ---------------------------
THREADS=16
HMM_OUT_MAG="~/analisis/hmmer/mags"
PHMMER_OUT_MAG="~/analisis/phmmer/mags"
HMM_OUT_REG="~/analisis/hmmer/reg"
PHMMER_OUT_REG="~/analisis/phmmer/reg"
HOM_ALIN="~/analisis/hom_alin"
TARGET_DIR="~/analisis/hom_hit"
POR_FAMILIA="~/analisis/hom_por_familia"


# ========================
# DATA CELANSING
# ========================

HOMOLOGY_HIT="~/analisis/hom_hit_"
DOMINIOS="~/analisis/dominios"
mkdir -p  "$HOMOLOGY_HIT" "$DOMINIOS"

# MOVE _hit.faa
find "$PHMMER_OUT_MAGS/" -type f -name '*_hits.faa' -exec cp {} "$HOMOLOGY_HIT" \;
find "$HMM_OUT_MAGS/"   -type f -name '*_hits.faa' -exec cp {} "$HOMOLOGY_HIT" \;
find "$PHMMER_OUT_REF/" -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_HIT" \;
find "$HMM_OUT_REF/"   -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_HIT" \;

# MAKE LIST OF NAMES
for i in "$HOMOLOGY_HIT"/*.faa; do
  basename "$i" .faa
done > "$HOMOLOGY_HIT/list.txt"

# MOVE _tab.txt
HOMOLOGY_TAB="~/analisis/hom_tab"

mkdir -p "$HOMOLOGY_TAB"
find "$PHMMER_OUT_MAGS/" -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_TAB" \;
find "$HMM_OUT_MAGS/"   -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_TAB" \;
find "$PHMMER_OUT_REF/" -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_TAB" \;
find "$HMM_OUT_REF/"   -type f -name '*_tab.txt' -exec cp {} "$HOMOLOGY_TAB" \;


# Filter homologs by e-value
# list with homologs data
find "$HOMOLOGY_TAB" -type f -name '*_tab.txt' -print0 \
  | xargs -0 grep -hv "^#" > "$HOMOLOGY_TAB"/all_hits_raw.txt

for f in "$HOMOLOGY_TAB"/*_tab.txt; do
    # Filters lines without comments and with e-value < 0.001, saves the file name
    awk -v file="$(basename "$f")" '!/^#/ && $5 > 0.001 {print file}' "$f" >> "$HOMOLOGY_TAB"/files_evalue_gt_0.001.txt
done


# Filter
TARGET_DIR="~/analisis/hom_hit"
mkdir -p "$TARGET_DIR"

while read -r line; do
    # Generate pattern of the corresponding .faa
    faa_pattern="${line/_res_tab.txt/_hits.faa}"


    # Search and copy
    find "$HOMOLOGY_HIT" -type f -name "$faa_pattern" -exec cp {} "$TARGET_DIR/" \;
done < "$HOMOLOGY_TAB/files_evalue_gt_0.001.txt"

# ========================
# ALIGN BY FAMILY
# ========================

# Extract family names (second part between "_" and before "_aligned_hits.faa")

FAMILIAS_MAG="~/analisis/familias_mag.txt"
FAMILIAS_REF="~/analisis/familias_ref.txt"
ID_FILE="~/analisis/list_no_fa_mags.txt"
ID_FILE_="~/analisis/list_no_fa_ref.txt"

sed -i 's/\r$//' "$ID_FILE"
sed -i 's/\r$//' "$ID_FILE_"

# Format correction
for file in "$TARGET_DIR"/*.faa; do
    [[ -f "$file" ]] || continue
    base="${file%.faa}"
    newbase=$(echo "$base" | tr ' ' '_' | tr -d '();' | tr '.' '_')
    newfile="${newbase}.faa"
    if [[ "$file" != "$newfile" ]]; then
        mv "$file" "$newfile"
        echo "Renamed: $file -> $newfile"
    fi
done


# Same wiht names
for f in "$ID_FILE" "$ID_FILE_"; do
    sed -i 's/^\s*//;s/\s*$//;s/\r//;s/[();]//g; s/\./_/g' "$f"
done


# Extract list of mechanisms in MAGs
while IFS= read -r id; do
  [[ -z "$id" ]] && continue
  for file in "$TARGET_DIR/${id}"_*.faa; do
    [[ -f "$file" ]] || continue
    base=$(basename "$file" .faa)
    familia=${base#${id}_}
    echo "$familia" >> "$FAMILIAS_MAG"
  done
done < "$ID_FILE"


# Extract list of mechanisms in Ref
while IFS= read -r id; do
  [[ -z "$id" ]] && continue
  for file in "$TARGET_DIR"/*"${id}"*.faa; do
    [[ -f "$file" ]] || continue
    base=$(basename "$file" .faa)
    familia=${base#${id}_}
    echo "$familia" >> "$FAMILIAS_REF"
  done
done < "$ID_FILE_"

POR_FAMILIA="~/analisis/hom_por_familia"
mkdir -p "$POR_FAMILIA"

# concatenate MAG families
while IFS= read -r fam; do
  [[ -z "$fam" ]] && continue
  out="$POR_FAMILIA/${fam}.faa"
  > "$out"
  for f in "$TARGET_DIR"/*"${fam}.faa"; do
    [[ -f "$f" ]] || continue
    # rename sequences by adding the file name as a prefix
    awk -v prefix="$(basename "$f" .faa)" '
      /^>/ {print ">" prefix "_" substr($0,2); next}
      {print}' "$f" >> "$out"
  done
done < "$FAMILIAS_MAG"


# concatenate REF families
while IFS= read -r fam; do
  [[ -z "$fam" ]] && continue
  out="$POR_FAMILIA/${fam}.faa"
  > "$out"
  for f in "$TARGET_DIR"/*"${fam}.faa"; do
    [[ -f "$f" ]] || continue
    # rename sequences by adding the file name as a prefix
    awk -v prefix="$(basename "$f" .faa)" '
      /^>/ {print ">" prefix "_" substr($0,2); next}
      {print}' "$f" >> "$out"
  done
done < "$FAMILIAS_REF"


HOM_ALIN="~/analisis/hom_alin"
mkdir -p "$HOM_ALIN"

# align by family
for f in "$POR_FAMILIA/"*.faa; do
  base=$(basename "$f" .faa)
  conda run -n align_env mafft --auto "$f" \
  > "$HOM_ALIN/${base}_aligned.faa"
done


# ================================
# HOMOLOGOUS PHYLOGENETIC TREES
# ================================

# cleaning with trimAL
HOM_ALIN_TRIM="~/analisis/hom_alin_trimmed"
mkdir -p "$HOM_ALIN_TRIM"

for aln in "$HOM_ALIN/"*.faa; do
  base=$(basename "$aln" _aligned.faa)

  conda run -n Hom_validation trimal \
    -in "$aln" \
    -out "$HOM_ALIN_TRIM/${base}_trimmed.faa" \
    -gt 0.7
done

# phylogenetic tree
HOM_ALIN_TREE="~/analisis/hom_tree"
mkdir -p "$HOM_ALIN_TREE"

for aln in "$HOM_ALIN_TRIM/"*.faa; do
  base=$(basename "$aln" .faa)
  conda run -n tree_env iqtree \
    -s "$aln" \
    -m LG+G+F \
    -nt $THREADS \
    -bb 1000 \
    -alrt 1000 \
    -pre "$HOM_ALIN_TREE/${base}"|| {
      echo "IQ-TREE fail with $aln, next..."
      continue
    }
done
```
