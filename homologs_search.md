# Bakta annotation: MAGS

```bash
#!/bin/bash

# CONFIGURATION
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# VARIABLES
ID_FILE="~/samples/id.txt"
BAKTA_OUT="~/binning/bakta"
BAKTA_DB="~/databases/bakta/db"
BINS_FILT_ALL="~/binning/filtered/all"
BINS_FILTER="~/binning/filtered/checkm2"

# MAKE DIRECTORIES
for dir in "$BAKTA_OUT" "$BINS_FILT_ALL" "$BINS_FILTER" "$BAKTA_DB"
do
    mkdir -p "$dir"
done

# ========================
# BAKTA ANNOTATION
# ========================

shopt -s nullglob

    while read -r ID; do
        for file in "$BINS_FILTER/$ID"/*.fa; do
            base=$(basename "$file")
            mv "$file" "$BINS_FILT_ALL/${ID}_$base"
        done  
    done < "$ID_FILE"

  # MAKE LIST.TXT
    for file in "$BINS_FILT_ALL"/*.fa; do
      basename "$file"
    done >> "$BINS_FILT_ALL/list.txt"

  # ANNOTATION
    while read -r line; do
      MAG_PATH="$BINS_FILT_ALL/$line"
      MAG_NAME=$(basename "$line" .fa)

      conda run -n bakta_env bakta \
        --db "$BAKTA_DB" \
        --prefix "$MAG_NAME" \
        --output "$BAKTA_OUT/$line" \
        --threads "$THREADS_PER_JOB" \
        "$MAG_PATH"
    done < "$BINS_FILT_ALL/list.txt"

  # MAKE LIST_NO_FA.TXT

    for id in "$BINS_FILT_ALL"/*.fa; do
      name=$(basename "$id" .fa)
      echo "$name"
    done > "$BINS_FILT_ALL/list_no_fa.txt"

```

# Bakta annotation: REF GENOMES
### In order to annotate the reference genomes, it is necessary to download and configure the card database.

```bash
# card
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
mkdir -p wildcard
tar -xjf wildcard_data.tar.bz2 -C wildcard
gunzip wildcard/*.gz


rgi card_annotation -i ~/databases/card_all/card.json > card_annotation.log 2>&1
rgi wildcard_annotation -i wildcard --card_json ~/databases/card_all/card.json -v 4.0.1 > wildcard_annotation.log 2>&1

rgi load \
  --card_json ~/databases/card_all/card.json \
  --debug --local \
  --card_annotation card_database_v4.0.1.fasta \
  --card_annotation_all_models card_database_v4.0.1_all.fasta \
  --wildcard_annotation wildcard_database_v4.0.1.fasta \
  --wildcard_annotation_all_models wildcard_database_v4.0.1_all.fasta \
  --wildcard_index ~/databases/card_all/wildcard/index-for-model-sequences.txt \
  --wildcard_version 4.0.1 \
  --amr_kmers ~/databases/card_all/wildcard/all_amr_61mers.txt \
  --kmer_database ~/databases/card_all/wildcard/61_kmer_db.json \
  --kmer_size 61
```
# Anotation

```bash
#!/bin/bash
set -euo pipefail

# CONFIGURATION
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# DIRECTORIES
SCRIPTS="~/scripts"
ID_FILE="~/binning/filtered/all/list_no_fa.txt"

BAKTA_OUT="~/binning/bakta"
BAKTA_OUT_REF="~/binning/bakta/ref"
BAKTA_DB="~/databases/bakta/db"

RGI_OUT="~/analisis/rgi"
CARD_FIL="~/databases/card"

# MAKE DIRECTORIES
for dir in "$SCRIPTS" "$BAKTA_OUT" "$BAKTA_DB" "$RGI_OUT" "$BAKTA_OUT_REF" \
"$CARD_FIL"; do
  mkdir -p "$dir"
done

# ========================
# EXTRACT CARD MECHANISM
# ========================

cd "$CARD_FIL"
tar -xvjf data ./protein_fasta_protein_homolog_model.fasta ./aro_index.tsv
mkdir "$CARD_FIL"/filt/
cd "$CARD_FIL"
python3 "$SCRIPTS"/extract_ARO.py

# Obtain bacterial genomes used in CARD
cd "$CARD_FIL"
conda run -n biopython python "$SCRIPTS"/get_card_genomes.py

# Filter single and complete species
awk -F '\t' 'NF && $2 ~ /^[A-Za-z]+ [A-Za-z]+$/ { print }' taxonomic_info.tsv | sort -t $'\t' -k2,2 -u > filtrado.tsv

# download species genomes
cd "$CARD_FIL"
cut -f2 filtrado.tsv > organisms.txt
mkdir -p "$CARD_FIL/card_sp"
cd "$CARD_FIL/card_sp"
# manual download

# MAKE LIST.TXT
for file in "$CARD_FIL"/card_sp/*.fna; do
  basename "$file"
done > "$CARD_FIL/card_sp/list_fa.txt"

# MAKE LIST_NO_FA.TXT
  for id in "$CARD_FIL"/card_sp/*.fna; do
    name=$(basename "$id" .fa)
    echo "$name"
  done > "$CARD_FIL/card_sp/list_no_fa.txt"

# ========================
# BAKTA ANNOTATION
# ========================

  while read -r line; do
    MAG_PATH="$CARD_FIL/card_sp/$line"
    MAG_NAME=$(basename "$line" .fa)

    conda run -n bakta_env bakta \
      --db "$BAKTA_DB" \
      --prefix "$MAG_NAME" \
      --output "$BAKTA_OUT_REF/$line" \
      --threads "$THREADS_PER_JOB" \
      "$MAG_PATH"
  done < "$CARD_FIL/card_sp/list_fa.txt" 
```

# Homolog search
```bash
#!/bin/bash
set -euo pipefail

# Configuracion
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# Rutas
SCRIPTS="~/scripts"
ID_FILE="~/binning/filtered/all/list_no_fa.txt"

BAKTA_OUT="~/binning/bakta"
BAKTA_OUT_REF="~/binning/bakta/ref"

RGI_OUT="~/analisis/rgi"
HMM_OUT_MAGS="~/analisis/hmmer/mags"
HMM_OUT_REF="~/analisis/hmmer/ref"
ALIGN_OUT="~/analisis/align"
TREE_OUT="~/analisis/tree"
CARD_FIL="~/databases/card"
PHMMER_OUT_MAGS="~/analisis/phmmer/mags"
PHMMER_OUT_REF="~/analisis/phmmer/ref"

HMM_B="~/analisis/hmmer/hmmbuild"


# Crear carpetas de salida
for dir in "$SCRIPTS" "$BAKTA_OUT" "$RGI_OUT" "$HMM_OUT_MAGS" "$ALIGN_OUT" "$BAKTA_OUT_REF" \
"$TREE_OUT" "$CARD_FIL" "$PHMMER_OUT_MAGS" "$HMM_B" "$HMM_OUT_REF" "$PHMMER_OUT_REF"; do
  mkdir -p "$dir"
done

mkdir -p "$ALIGN_OUT/multi"
mkdir -p "$ALIGN_OUT/single"

# =========================
# CARD MECHANISMS ALIGNMENT
# =========================

for file in "$CARD_FIL"/filt/*.fasta; do
    base_name=$(basename "$file" .fasta)
    count=$(grep -c "^>" "$file")


    if [ "$count" -ge 2 ]; then

        conda run -n align_env mafft \
        --auto \
        --thread "$THREADS_PER_JOB" \
        "$file" > "$ALIGN_OUT/multi/${base_name}_aligned.fasta"
 
    else

        short_name=$(echo "$base_name" | tr ';' '_' | cut -c1-50)

        mv "$file" "$ALIGN_OUT/single/${short_name}_aligned.fasta"
        echo "Movido $file como ${short_name}_aligned.fasta"
        
    fi
    
done

# ============================
# HMM BUILD OF CARD MECHANISMS
# ============================

mkdir -p "$HMM_B/multi"
for file in "$ALIGN_OUT"/multi/*_aligned.fasta; do
    base_name=$(basename "$file" .fasta)

    conda run -n hmmer_env hmmbuild \
        --cpu "$THREADS_PER_JOB" \
        "$HMM_B/multi/${base_name}.hmm" \
        "$file"
done

# HOMOLOGY SEARCH REF
# ========================
# HMM SEARCH
# ========================

run_hmmsearch() {
    file="$1"
    line="$2"

    base_name=$(basename "$file" .hmm)

    mkdir -p "$HMM_OUT_REF/$line"

    conda run -n hmmer_env hmmsearch \
        -o "$HMM_OUT_REF/$line/${base_name}.txt" \
        --tblout "$HMM_OUT_REF/$line/${line}_${base_name}_res_tab.txt" \
        --cpu "$THREADS_PER_JOB" \
        "$file" \
        "$BAKTA_OUT_REF"/"$line".fna/"$line".faa

    # Filter e-value < 1e-5
    awk '!/^#/ && $5 < 1e-5 {print $1}' "$HMM_OUT_REF/$line/${line}_${base_name}_res_tab.txt" > \
        "$HMM_OUT_REF/$line/${line}_${base_name}_hits.txt"

    # Extract matches
    conda run -n hmmer_env seqtk subseq \
        "$BAKTA_OUT_REF"/"$line".fna/"$line".faa "$HMM_OUT_REF/$line/${line}_${base_name}_hits.txt" > \
        "$HMM_OUT_REF/$line/${line}_${base_name}_hits.faa"
}

# Export variables and function
export HMM_B HMM_OUT_REF BAKTA_OUT_REF ID_FILE THREADS_PER_JOB NUM_PARALLEL_JOBS
export -f run_hmmsearch

# Execute in paralell
parallel --jobs "$NUM_PARALLEL_JOBS" run_hmmsearch ::: "$HMM_B"/multi/*.hmm ::: $(cat "$CARD_FIL/card_sp/list_no_fa.txt")


# ========================
# PHMMER SEARCH
# ========================

run_phmmer() {
    file="$1"
    line="$2"

    base_name=$(basename "$file" .fasta)

    mkdir -p "$PHMMER_OUT_REF/$line"

    conda run -n hmmer_env phmmer \
        -o "$PHMMER_OUT_REF/$line/${base_name}.txt" \
        --tblout "$PHMMER_OUT_REF/$line/${line}_${base_name}_res_tab.txt" \
        --cpu "$THREADS_PER_JOB" \
        "$BAKTA_OUT_REF"/"$line".fna/"$line".faa \
        "$file"

    # FILTER e-VALUE < 1e-5 AND CLEAN ID
    awk '!/^#/ && $5 < 1e-5 { print $3 }' \
        "$PHMMER_OUT_REF/$line/${line}_${base_name}_res_tab.txt" > \
        "$PHMMER_OUT_REF/$line/${line}_${base_name}_hits.txt"            

    # EXTRACT MATCHES
    conda run -n hmmer_env seqtk subseq \
        "$BAKTA_OUT_REF"/"$line".fna/"$line".faa "$PHMMER_OUT_REF/$line/${line}_${base_name}_hits.txt" > \
        "$PHMMER_OUT_REF/$line/${line}_${base_name}_hits.faa"
}

# Export environment and function
export PHMMER_OUT_REF BAKTA_OUT_REF ALIGN_OUT ID_FILE THREADS_PER_JOB
export -f run_phmmer

# Execute in paralell
parallel --jobs "$NUM_PARALLEL_JOBS" run_phmmer ::: "$ALIGN_OUT"/single/*.fasta ::: $(cat "$CARD_FIL/card_sp/list_no_fa.txt")


# HOMOLOGY SEARCH MAGS
# ========================
# HMM SEARCH
# ========================

run_hmmsearch() {
    file="$1"
    line="$2"

    base_name=$(basename "$file" .hmm)

    mkdir -p "$HMM_OUT_MAGS/$line"

    conda run -n hmmer_env hmmsearch \
        -o "$HMM_OUT_MAGS/$line/${base_name}.txt" \
        --tblout "$HMM_OUT_MAGS/$line/${line}_${base_name}_res_tab.txt" \
        --cpu "$THREADS_PER_JOB" \
        "$file" \
        "$BAKTA_OUT/$line.fa/$line.faa"

    # Filter e-value < 1e-5
    awk '!/^#/ && $5 < 1e-5 {print $1}' "$HMM_OUT_MAGS/$line/${line}_${base_name}_res_tab.txt" > \
        "$HMM_OUT_MAGS/$line/${line}_${base_name}_hits.txt"

    # Extract matches
    conda run -n hmmer_env seqtk subseq \
        "$BAKTA_OUT/$line.fa/$line.faa" "$HMM_OUT_MAGS/$line/${line}_${base_name}_hits.txt" > \
        "$HMM_OUT_MAGS/$line/${line}_${base_name}_hits.faa"
}

# Export variables and funcition
export HMM_B HMM_OUT_MAGS BAKTA_OUT ID_FILE THREADS_PER_JOB NUM_PARALLEL_JOBS
export -f run_hmmsearch

# Execute in paralell
parallel --jobs "$NUM_PARALLEL_JOBS" run_hmmsearch ::: "$HMM_B"/multi/*.hmm ::: $(cat "$ID_FILE")



# ========================
# PHMMER SEARCH
# ========================

run_phmmer() {
    file="$1"
    line="$2"

    base_name=$(basename "$file" .fasta)

    mkdir -p "$PHMMER_OUT_MAGS/$line"

    conda run -n hmmer_env phmmer \
        -o "$PHMMER_OUT_MAGS/$line/${base_name}.txt" \
        --tblout "$PHMMER_OUT_MAGS/$line/${line}_${base_name}_res_tab.txt" \
        --cpu "$THREADS_PER_JOB" \
        "$BAKTA_OUT/$line.fa/$line.faa" \
        "$file"

    # FILTER e-VALUE < 1e-5 AND CLEAN ID
    awk '!/^#/ && $5 < 1e-5 { print $3 }' \
        "$PHMMER_OUT_MAGS/$line/${line}_${base_name}_res_tab.txt" > \
        "$PHMMER_OUT_MAGS/$line/${line}_${base_name}_hits.txt"            

    # EXTRACT MATCHES
    conda run -n hmmer_env seqtk subseq \
        "$BAKTA_OUT/$line.fa/$line.faa" "$PHMMER_OUT_MAGS/$line/${line}_${base_name}_hits.txt" > \
        "$PHMMER_OUT_MAGS/$line/${line}_${base_name}_hits.faa"
}

# Export enviroment and function
export PHMMER_OUT_MAGS BAKTA_OUT ALIGN_OUT ID_FILE THREADS_PER_JOB
export -f run_phmmer

# Execute in parallel
parallel --jobs "$NUM_PARALLEL_JOBS" run_phmmer ::: "$ALIGN_OUT"/single/*.fasta ::: $(cat "$ID_FILE")
```
