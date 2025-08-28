# PHYLOGENETIC TREES OF SPECIES

```bash

#!/bin/bash
set -euo pipefail

# ---------------------------
# CONFIGURATION
# ---------------------------
FAAS="~/analisis/tree/bins"   
BUSCO_LINEAGE="~/databases/busco/bacteria_odb12"
BUSCO_OUT="~/analisis/tree/busco"
THREADS=20
ALIGN_DIR="~/analisis/tree/alignments"
CONCAT_ALIGNMENT="~/analisis/tree"
TREE_OUT="~/analisis/tree/phylo_tree"
GENOMES_DIR="~/analisis/tree/genomes"
PROKKA_OUT="~/analisis/tree/prokka_out"
MARKER_DIR="$BUSCO_OUT/markers_by_gene"
TRIM_DIR="~/analisis/tree/alignments_trimmed"

mkdir -p "$TRIM_DIR"
mkdir -p "$MARKER_DIR"
mkdir -p "$FAAS"
mkdir -p "$ALIGN_DIR"
mkdir -p "$TREE_OUT"
mkdir -p "$PROKKA_OUT"

# Out-group
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/465/GCF_000011465.1_ASM1146v1/GCF_000011465.1_ASM1146v1_genomic.fna.gz -O Aquifex_aeolicus.fna.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/545/GCF_000008545.1_ASM854v1/GCF_000008545.1_ASM854v1_genomic.fna.gz -O Thermotoga_maritima.fna.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/585/GCF_000008585.1_ASM858v1/GCF_000008585.1_ASM858v1_genomic.fna.gz -O Deinococcus_radiodurans.fna.gz


# ---------------------------
# # Genomic annotation PROKKA
# ---------------------------


for GENOME in "$GENOMES_DIR"/*.fna; do
    BASENAME=$(basename "$GENOME" .fna)
    prokka --outdir "$PROKKA_OUT/$BASENAME" \
           --prefix "$BASENAME" \
           --cpus 10 \
           "$GENOME"
done


# ---------------------------
# CLEAN NAMES
# ---------------------------


find ~/analisis/tree/prokka_out -type f -name "*.faa" | while read -r f; do
  name=$(basename "$f" .faa)
  clean_name=$(echo "$name" | tr '.' '_')
  cp "$f" ~/analisis/tree/bins/"$clean_name".faa
done



# ---------------------------
# RUN BUSCO
# ---------------------------
for f in "$FAAS"/*.faa; do
    name=$(basename "$f" .faa)
    conda run -n eval_bins_env_busco busco -i "$f" -l "$BUSCO_LINEAGE" -m prot -o "busco_$name" -c "$THREADS"
done


# ---------------------------
# REORDER BY BOOKMARK
# ---------------------------

for d in $BUSCO_OUT/busco_*/run_bacteria_odb12/busco_sequences/single_copy_busco_sequences/; do
    GENOME=$(basename $(dirname $(dirname $(dirname $d))))  # extact busco_<genome>
    
    for marker in "$d"/*.faa; do
        gene=$(basename "$marker" .faa)
        mkdir -p "$MARKER_DIR/$gene"


        # Rename the sequence within the file to the genome name
        sed "s/>.*/>$GENOME/" "$marker" > "$MARKER_DIR/$gene/$GENOME.faa"
    done
done

# Concatenate all bookmarks into a single file per bookmark
for gene_dir in "$MARKER_DIR"/*; do
    gene=$(basename "$gene_dir")
    cat "$gene_dir"/*.faa > "$MARKER_DIR/${gene}_unido.faa"
done



# -----------------------------
# ALIGN EACH MARKER WITH MAFFT
# -----------------------------
for gene_file in "$MARKER_DIR"/*_unido.faa; do
    gene=$(basename "$gene_file" .faa)
    conda run -n eval_bins_env_busco mafft --maxiterate 1000 --localpair --thread 10 --auto "$gene_file" > "$ALIGN_DIR/$gene.aln"
done


# ---------------------------------
# FILTER ALIGMENTS WITH trimAl
# ---------------------------------

for aln in "$ALIGN_DIR"/*.aln; do
    base=$(basename "$aln" .aln)
    conda run -n Hom_validation trimal -in "$aln" -out "$TRIM_DIR/$base.trim.aln" -gt 0.7
done


# ---------------------------
# CONCATENATE ALIGMENTS
# ---------------------------
# Using AMAS
conda run -n amas AMAS.py concat -i "$TRIM_DIR"/*.aln -f fasta -d aa -p "$CONCAT_ALIGNMENT"


# ---------------------------
# BUILD A TREE WITH IQ-TREE
# ---------------------------

# base tree
conda run -n tree_env iqtree \
  -s "/home/omar/Proyecto/analisis/tree/busco/concatenated.out" \
  -m LG+G+F \
  -nt "$THREADS" \
  -pre "$TREE_OUT/phylo_tree"

# support
conda run -n tree_env iqtree \
  -s "/home/omar/Proyecto/analisis/tree/busco/concatenated.out" \
  -m LG+G+F \
  -nt "$THREADS" \
  -bb 1000 \
  -alrt 1000 \
  -t "$TREE_OUT/phylo_tree.treefile" \
  -redo \
  -pre "$TREE_OUT/phylo_tree_ufboot"
```
