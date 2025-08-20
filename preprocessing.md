# Raw reads preprocessing


```bash

# ---------------------------
# First QC
# ---------------------------

# .fastq in samples directory
conda run -n metabiome-preprocessing fastqc *.fastq

# ---------------------------
# Adapters filter
# ---------------------------

for file in ~/samples/*_1.fastq; do
  basename "$file" _1.fastq
done > ~/samples/id.txt

while read line; do
    conda run -n metabiome-preprocessing trim_galore --cores 20 --phred33 --length 100 --stringency 3 --paired \
    -o ~/trimgalore \
    ~/samples/"$line"_1.fastq ~/samples/"$line"_2.fastq
done < ~/samples/id.txt

# ---------------------------
# Remove host sequence
# ---------------------------
# PRJNA1269778: Download sequences of Homo sapiens, Rinolophus ferromequinum, Miniopterus schrebersi, Myotis myotis.
# PRJNA954561: Download sequences of Homo sapiens, Rhinolophus sinicus, Myotis pilosus, Miniopterus pusillus, Murina feae,
# Rousettus leschenaultii, Lyroderma lyra, Rhinolophus affinis, Hipposideros armiger, Miniopterus fuliginosus,
# Eonycteris spelaea, Hipposideros larvatus, Pipistrellus abramus, Eptesicus fuscus, Myotis davidii.

# Join genomes
cat ~/bowtie2/ind_bow/*.fnas  > Mixed.fasta

# Remove host
while read sample_id; do
  bowtie2 \
    -x ~/bowtie2/ind_bow/Mix \
    -1 ~/trimgalore/${sample_id}_1_val_1.fq \
    -2 ~/trimgalore/${sample_id}_2_val_2.fq \
    -q \
    -p 20 \
    --un-conc-gz ~/bowtie2/no_pareados/${sample_id}_unmapped.fq.gz \
    -S /dev/null
  echo "Processed sample: $sample_id"
done < ~/samples/id.txt


# ---------------------------
# Second QC
# ---------------------------
conda run -n metabiome-preprocessing fastqc *.fastq
```
