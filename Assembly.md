# Assembly

#### In order to proceed with the assembly, it is necessary to download the busco, checkm2 and bakta databases.
```bash
# checkm2
checkm2 database --download --path ~/databases/checkm_data
export CHECKM2DB="~/databases/checkm_data"

# busco
wget https://busco-data.ezlab.org/v5/data/lineages/fungi_odb12.2025-07-01.tar.gz
tar -xvzf fungi_odb12.2025-07-01.tar.gz

# BUSCO_LINEAGE="eukaryota_odb12"
# Check disponibility: busco --list-datasets

| BUSCO             | Taxonomy  |
| ----------------- | ----------|
| `eukaryota_odb12` | Eucariota |
| `bacteria_odb12`  | Bacteria  |
| `archaea_odb12`   | Arquea    |
| `viruses_odb12`   | Virus     |

# bakta
bakta_db download --output <output-path> --type [light|full]
```

```bash
#!/bin/bash

# CONFIGURATION
NUM_PARALLEL_JOBS=1
THREADS_PER_JOB=20

# DIRECTORIES
SCRIPTS="~/scripts"
UNMAPPED_DIR="~/bowtie2/no_pareados"
UNAS_DIR="~/bowtie2/unassembled"
MAPPED_MAG_DIR="~/bowtie2/mags_mapped"
DEDUPLI_DIR="~/fastp"
OUT_DIR_BASE1="~/megahit"
OUT_DIR_BASE="~/metaspades"
ID_FILE="~/samples/id.txt"
IND_BOW="~/ind_bow/megahit_salamanca"
MAGS="~/mags"
IND_MAGS="~/ind_bow/mags"
BIN_DATA="~/binning"
MET_BIN="~/binning/metabat/bin"
MAX_BIN="~/binning/maxbin/bin"
CONCOCT_BIN="~/binning/concoct"
DASTOOL_OUT="~/binning/dastool"
CHECKM2_OUT="~/binning/checkm2"
BUSCO_OUT="~/binning/busco"
CHECKM2_DB="~/databases/checkm_data/CheckM2_database"
BUSCO_DB="~/databases/busco_database"
BAKTA_DB="~/databases/bakta/db"
BINS_FILTER="~/binning/filtered/checkm2"
BINS_FILT_BU="~/binning/filtered/busco"
BAKTA_OUT="~/binning/bakta"
BINS_FILT_ALL="~/binning/filtered/all"

# MAKE DIRECTORIES
for dir in "$OUT_DIR_BASE1" "$OUT_DIR_BASE" "$IND_BOW" "$UNAS_DIR" "$BUSCO_DB" \
"$MAGS" "$IND_MAGS" "$MAPPED_MAG_DIR" "$MET_BIN" "$MAX_BIN" "$BIN_DATA" "$BINS_FILTER" \
"$CONCOCT_BIN" "$DASTOOL_OUT" "$SCRIPTS" "$CHECKM2_OUT" "$BUSCO_OUT" "$CHECKM2_DB" \
"$BINS_FILT_BU" "$BAKTA_DB" "$BAKTA_OUT" "$BINS_FILT_ALL"
do
    mkdir -p "$dir"
done


# EXPORT TO PARALLEL
export DEDUPLI_DIR OUT_DIR_BASE1 THREADS_PER_JOB IND_BOW UNMAPPED_DIR DASTOOL_OUT SCRIPTS \
UNAS_DIR OUT_DIR_BASE MAGS IND_MAGS MAPPED_MAG_DIR BIN_DATA MET_BIN MAX_BIN CONCOCT_BIN \
CHECKM2_OUT BUSCO_OUT CHECKM2_DB BUSCO_DB BINS_FILTER BINS_FILT_BU

# EXECUTE IN PARALLEL
parallel -j "$NUM_PARALLEL_JOBS" --bar '
  set -euo pipefail

  # MEGAHIT
    conda run -n ensamblaje_env megahit \
        -1 "$DEDUPLI_DIR/{}_1.fastq.gz" \
        -2 "$DEDUPLI_DIR/{}_2.fastq.gz" \
        --out-dir "$OUT_DIR_BASE1/{}" \
        --num-cpu-threads $THREADS_PER_JOB
    echo "Processed by megahit sample: {}"

  # FILTER MEGAHIT >1500
    conda run -n ensamblaje_env seqkit seq \
      -m 1500 \
      -g "$OUT_DIR_BASE1/{}/final.contigs.fa" \
      > "$OUT_DIR_BASE1/{}/contigs_mayores_1500.fa"
    echo "Filtered >1500 sample: {}"

  # INDEX AND MAP BOWTIE2
      conda run -n metabiome-preprocessing bowtie2-build \
      "$OUT_DIR_BASE1/{}/contigs_mayores_1500.fa" \
      "$IND_BOW/{}"
    echo "Indexed sample: {}"

    conda run -n metabiome-preprocessing bowtie2 \
      -x "$IND_BOW/{}" \
      -1 "$UNMAPPED_DIR/{}_unmapped_1.fastq.gz" \
      -2 "$UNMAPPED_DIR/{}_unmapped_2.fastq.gz" \
      -q \
      -p $THREADS_PER_JOB \
      --un-conc-gz "$UNAS_DIR/tmp_{}_unmapped" \
      -S /dev/null
    echo "Processed by bowtie2 sample: {}"

  # RENAME UNMAPPED
    mv "$UNAS_DIR/tmp_{}_unmapped.1" \
      "$UNAS_DIR/{}_unmapped_1.fastq.gz"

    mv "$UNAS_DIR/tmp_{}_unmapped.2" \
      "$UNAS_DIR/{}_unmapped_2.fastq.gz"
    echo "Renamed sample: {}"


# CORRECTION WITH FASTP FOR SAMPLES 90, 91 Y 92
conda run -n ensamblaje_env fastp \
  -i "$UNAS_DIR/{}_unmapped_1.fastq.gz" \
  -I "$UNAS_DIR/{}_unmapped_2.fastq.gz" \
  -o "$UNAS_DIR/{}_clean_1.fastq.gz" \
  -O "$UNAS_DIR/{}_clean_2.fastq.gz" \
  --detect_adapter_for_pe \
  --trim_front1 10 --trim_front2 10 \
  --thread $THREADS_PER_JOB


  #DIGITAL NORMALIZATION
    conda run -n ensamblaje_env bbnorm.sh \
      in="$UNAS_DIR/{}_clean_1.fastq.gz" \
      in2="$UNAS_DIR/{}_clean_2.fastq.gz" \
      out="$UNAS_DIR/{}_norm_1.fastq" \
      out2="$UNAS_DIR/{}_norm_2.fastq" \
      target=40 min=5 threads=$THREADS_PER_JOB passes=2
    echo "Normalized sample: {}"

  # METASPADES
    conda run -n ensamblaje_env metaspades.py \
      -k 21,33 \
      -1 "$UNAS_DIR/{}_norm_1.fastq" \
      -2 "$UNAS_DIR/{}_norm_2.fastq" \
      -o "$OUT_DIR_BASE/{}/" \
      -t $THREADS_PER_JOB \
      -m 44
    echo "Processed by metaspades sample: {}"
    
  # FILTER METASPADES >1500
    conda run -n ensamblaje_env seqkit seq \
      -m 1500 \
      -g "$OUT_DIR_BASE/{}/contigs.fasta" \
      > "$OUT_DIR_BASE/{}/contigs_mayores_1500.fa"
    echo "Filtered >1500 metaspades sample: {}"

  # RENAME PREFIXES
    sed "s/^>/>megahit_{}_/" \
      "$OUT_DIR_BASE1/{}/contigs_mayores_1500.fa" \
      > "$OUT_DIR_BASE1/{}/{}_megahit_renamed.fa"
    echo "Renamed megahit sample: {}"

    sed "s/^>/>metaspades_{}_/" \
      "$OUT_DIR_BASE/{}/contigs_mayores_1500.fa" \
      > "$OUT_DIR_BASE/{}/{}_metaspades_renamed.fa"
    echo "Renamed metaspades sample: {}"

  # JOIN MAGS
    mkdir -p "$MAGS/{}"

    cat "$OUT_DIR_BASE1/{}/{}_megahit_renamed.fa" \
      "$OUT_DIR_BASE/{}/{}_metaspades_renamed.fa" \
      > "$MAGS/{}/{}_may_1500_unido.fa"
    echo "Unite megahit_metaspades completed: {}"

  # MAGS INDEX AND MAPED WITH BOWTIE2
    conda run -n metabiome-preprocessing bowtie2-build\
      "$MAGS/{}/{}_may_1500_unido.fa" \
      "$IND_MAGS/{}"
    echo "Indexed sample: {}"

    conda run -n metabiome-preprocessing bowtie2 \
      -x "$IND_MAGS/{}" \
      -1 "$UNMAPPED_DIR/{}_unmapped_1.fastq.gz" \
      -2 "$UNMAPPED_DIR/{}_unmapped_2.fastq.gz" \
      -q \
      -p $THREADS_PER_JOB \
      -S "$MAPPED_MAG_DIR/{}_mapped.sam"
    echo "Processed by bowtie2 sample: {}"

  # CONVERT TO BAM
      conda run -n ensamblaje_env samtools \
      view -bS "$MAPPED_MAG_DIR/{}_mapped.sam" -o "$MAPPED_MAG_DIR/{}_mapped.bam"
    echo "converted a bam: {}"

  # SORT
    conda run -n ensamblaje_env samtools \
      sort "$MAPPED_MAG_DIR/{}_mapped.bam" \
      -o "$MAPPED_MAG_DIR/{}_sorted.bam"
    echo "Sortered: {}"
    
  # INDEX
    conda run -n ensamblaje_env samtools \
      index "$MAPPED_MAG_DIR/{}_sorted.bam"

  # CALCULATE DEPTH OF COVERAGE WITH JGI
    conda run -n ensamblaje_env jgi_summarize_bam_contig_depths \
      --outputDepth "$BIN_DATA/{}_depth.txt" \
      "$MAPPED_MAG_DIR/{}_sorted.bam"

  # METABAT2
    conda run -n ensamblaje_env metabat2 \
      -i "$MAGS/{}/{}_may_1500_unido.fa" \
      -a "$BIN_DATA/{}_depth.txt" \
      -o "$MET_BIN/{}_bin"

  # CALCULATE ABUNDANCE
    conda run -n ensamblaje_env jgi_summarize_bam_contig_depths \
      --outputDepth "$BIN_DATA/{}_abund.abund" \
      "$MAPPED_MAG_DIR/{}_sorted.bam"

  # MAXBIN
    conda run -n ensamblaje_env run_MaxBin.pl \
      -contig "$MAGS/{}/{}_may_1500_unido.fa" \
      -abund "$BIN_DATA/{}_abund.abund" \
      -out "$MAX_BIN/{}_bin"


  # FRAGMENT CONTIGS in 10 kb
    conda run -n ensamblaje_env cut_up_fasta.py \
      "$MAGS/{}/{}_may_1500_unido.fa" \
      -c 10000 \
      -o 0 \
      --merge_last \
      -b "$MAGS/{}/{}_contigs_cut10k.bed" > "$MAGS/{}/{}_contigs_cut10k.fa"

 
  # INDEX AND MAP
    conda run -n metabiome-preprocessing bowtie2-build \
      "$MAGS/{}/{}_may_1500_unido.fa" "$IND_MAGS/{}_contigs_cut"

    conda run -n metabiome-preprocessing bowtie2 \
      -x "$IND_MAGS/{}_contigs_cut" \
      -1 "$UNMAPPED_DIR/{}_unmapped_1.fastq.gz" \
      -2 "$UNMAPPED_DIR/{}_unmapped_2.fastq.gz" \
      -S "$MAPPED_MAG_DIR/{}_mapped_cut.sam" \
      -q \
      -p $THREADS_PER_JOB
    echo "Processed by bowtie2 sample: {}"

    conda run -n ensamblaje_env samtools \
    view -bS "$MAPPED_MAG_DIR/{}_mapped_cut.sam" -o "$MAPPED_MAG_DIR/{}_mapped_cut.bam"

    conda run -n ensamblaje_env samtools \
    sort "$MAPPED_MAG_DIR/{}_mapped_cut.bam" -o "$MAPPED_MAG_DIR/{}_mapped_cut_sorted.bam"

    conda run -n ensamblaje_env samtools \
    index "$MAPPED_MAG_DIR/{}_mapped_cut_sorted.bam"

  # COVERAGE BY POSITION
    conda run -n ensamblaje_env concoct_coverage_table.py \
      "$MAGS/{}/{}_contigs_cut10k.bed" \
      "$MAPPED_MAG_DIR/{}_mapped_cut_sorted.bam" \
      > "$BIN_DATA/{}_coverage_per_base.tsv"
    echo "Cobertura por posicion sample: {}"

  # CONCOCT
    conda run -n ensamblaje_env concoct \
      --composition_file "$MAGS/{}/{}_contigs_cut10k.fa" \
      --coverage_file "$BIN_DATA/{}_coverage_per_base.tsv" \
      -b "$CONCOCT_BIN/{}_bin"
    echo "Processed by CONCOCT sample: {}"

  # CONVERT TO .scaffolds2bin.tsv

    conda run -n ensamblaje_env "$SCRIPTS"/metabat2_to_tsv.sh "$MET_BIN"/{}_bin*.fa > "$MET_BIN/{}_metabat2.scaffolds2bin.tsv"
    echo "scaffolds2bin for metabat sample: {}"

    conda run -n ensamblaje_env "$SCRIPTS"/maxbin2_to_tsv.sh "$MAX_BIN"/{}_bin*.fasta > "$MAX_BIN/{}_maxbin2.scaffolds2bin.tsv"
    echo "scaffolds2bin for maxbin sample: {}"

  # MERGE OF CONCOCT FOR DAS_TOOL
    conda run -n ensamblaje_env merge_cutup_clustering.py \
      "$CONCOCT_BIN"/{}_bin_clustering_gt1000.csv \
      > "$CONCOCT_BIN"/{}_merged_clustering.tsv

    sed '\''1d; s/,/\t/g'\'' "$CONCOCT_BIN"/{}_merged_clustering.tsv > "$CONCOCT_BIN"/{}_concoct.scaffolds2bin.tsv
    echo "scaffolds2bin para Concoct sample: {}"


  # DAS_Tool
    mkdir -p "$DASTOOL_OUT/{}_bins"

    conda run -n ensamblaje_env DAS_Tool \
      -i "$MET_BIN/{}_metabat2.scaffolds2bin.tsv","$MAX_BIN/{}_maxbin2.scaffolds2bin.tsv","$CONCOCT_BIN/{}_concoct.scaffolds2bin.tsv" \
      -l metabat2,maxbin2,concoct \
      -c "$MAGS/{}/{}_may_1500_unido.fa" \
      -o "$DASTOOL_OUT/{}_bins/{}_dastool" \
      --search_engine diamond \
      --write_bins \
      --threads $THREADS_PER_JOB
    echo "Dastool for sample: {}"

  # BINS EVALUATION
    echo "Running CheckM2"
    conda run -n eval_bins_env_checkm checkm2 \
      predict \
      -i "$DASTOOL_OUT/{}_bins/{}_dastool_DASTool_bins" \
      -x fa \
      -o "$CHECKM2_OUT/{}_checkm2" \
      --database_path "$CHECKM2_DB/uniref100.KO.1.dmnd" \
      --threads $THREADS_PER_JOB

  # {} = ID
  # $THREADS_PER_JOB
  # $DASTOOL_OUT = rute to bins of DASTool
  # $BUSCO_OUT = rute to BUSCO results
  # fungi_odb12 = data base of BUSCO lineage

    conda run -n eval_bins_env_busco bash '"$SCRIPTS"'/run_busco.sh {} \
    '"$THREADS_PER_JOB"' '"$DASTOOL_OUT"' '"$BUSCO_OUT"' '"$BUSCO_DB/fungi_odb12"'

  # FILTER BINS CHECKM2
    awk -v min_comp=90 -v max_cont=5 -F'\''\t'\'' '\'' 
      NR>1 && $2 >= 90 && $3 <= 5 { print $1 }
    '\'' "$CHECKM2_OUT/{}_checkm2/quality_report.tsv" > "$CHECKM2_OUT/{}_checkm2/{}_bins_filtrados.txt"

  # COPY FILTERED
    mkdir -p "$BINS_FILTER/{}"

    while read -r bin_name; do
      cp "$DASTOOL_OUT/{}_bins/{}_dastool_DASTool_bins/${bin_name}.fa" "$BINS_FILTER/{}/" || echo "No $bin_name"
    done < "$CHECKM2_OUT/{}_checkm2/{}_bins_filtrados.txt"

  # FILTER WITH BUSCO
    conda run -n eval_bins_env_busco bash '"$SCRIPTS"'/filter_busco.sh \
        '"${BUSCO_OUT}/${ID}"' '"${BINS_FILT_BU}/{}_busco_filtrados.txt"'


' :::: "$ID_FILE"
```
