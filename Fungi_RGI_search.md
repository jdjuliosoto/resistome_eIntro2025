# Fungi search

### In order to search for bins corresponding to mushrooms, it is necessary to download the kraken2 database.
```bash
mkdir -p ~/databases/kraken_fungi_db
kraken2-build --download-taxonomy --db ~/databases/kraken_fungi_db
kraken2-build --download-library fungi --db ~/databases/kraken_fungi_db
kraken2-build --build --db ~/databases/kraken_fungi_db


# Search fungi
while read line; do
  conda run -n fungi_env kraken2 \
    --db ~/databases/kraken_fungi_db \
    --paired ~/bowtie2/no_pareados/${line}_unmapped_1.fastq.gz \
    ~/bowtie2/no_pareados/${line}_unmapped_2.fastq.gz \
    --threads 20 \
    --report ~/Kraken_fungi/${line}_kraken2_report.txt \
    --classified-out ~/Kraken_fungi/classified/${line}_classified#_fungi.fastq \
    --unclassified-out ~/Kraken_fungi/unclassified/${line}_unclassified#_nonfungi.fastq \
    --output ~/Kraken_fungi/${line}_kraken2_output.txt
done < "~/samples/id.txt"

```
# RGI search
```bash
#!/bin/bash
set -euo pipefail

# Configuracion
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# Rutas
ID_FILE="~/binning/filtered/all/list_no_fa.txt"
ID_ORIG="~/samples/id.txt"
BAKTA_OUT="~/binning/bakta"
RGI_OUT="~/analisis/rgi"
MAG="~/mags"
RAW="~/bowtie2/no_pareados"

# Crear carpetas de salida
for dir in "$BAKTA_OUT" "$RGI_OUT" "$MAG" "$RAW" ; do
  mkdir -p "$dir"
done

# Export to parallel
export ID_FILE BAKTA_OUT RGI_OUT THREADS_PER_JOB NUM_PARALLEL_JOBS MAG RAW THREADS_PER_JOB

# Run in parallel for BINS
parallel -j "$NUM_PARALLEL_JOBS" --bar "
  set -euo pipefail

  SAMPLE_ID={}
  FAA_FILE=\"$BAKTA_OUT/\$SAMPLE_ID/\${SAMPLE_ID}.faa\"

  mkdir -p \"$RGI_OUT/bins\"

  # RESISTANCE CARD (RGI)
  conda run -n rgi_env rgi main \
    --input_sequence \"\$FAA_FILE\" \
    --output_file \"$RGI_OUT/bins/\${SAMPLE_ID}_rgi\" \
    --input_type protein \
    --num_threads $THREADS_PER_JOB \
    --clean
" :::: "$ID_FILE"


# Run in parallel for MAGs
parallel -j "$NUM_PARALLEL_JOBS" --bar "
  set -euo pipefail

  SAMPLE_ID={}
  MAG_FILE=\"$MAG/\$SAMPLE_ID/\${SAMPLE_ID}_may_1500_unido.fa\"

  mkdir -p \"$RGI_OUT/mags\"

  conda run -n rgi_env rgi main \
    --input_sequence \"\$MAG_FILE\" \
    --output_file \"$RGI_OUT/mags/\${SAMPLE_ID}_rgi\" \
    --input_type contig \
    --num_threads $THREADS_PER_JOB \
    --clean
" :::: "$ID_ORIG"


# Run RAWREADS in parallel
parallel -j "$NUM_PARALLEL_JOBS" --bar "
  set -euo pipefail

  export RGI_DATABASE=~/databases/card/localDB

  SAMPLE_ID={}
  RAW_FILE_1=\"$RAW/\${SAMPLE_ID}_unmapped_1.fastq.gz\"
  RAW_FILE_2=\"$RAW/\${SAMPLE_ID}_unmapped_2.fastq.gz\"

  mkdir -p \"$RGI_OUT/raws\"

  conda run -n rgi_env rgi bwt \
    --read_one \"\$RAW_FILE_1\" \
    --read_two \"\$RAW_FILE_2\" \
    --output_file \"$RGI_OUT/raws/\${SAMPLE_ID}_rgi\" \
    --local \
    --threads $THREADS_PER_JOB \
    --clean \
    --aligner bowtie2
" :::: "$ID_ORIG"
```
