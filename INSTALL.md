# Installation Instructions

## **Hardware and Software Setup**
All analyses were performed using a workstation equipped with:

* Processor : Intel(R) Core(TM) i7-14700 CPU @ 5.4 GHz (20 cores, 28 threads) 
* Main Drive (OS + apps):  
    - Type: NVMe SSD  
    - Capacity: 512 GB  
    - Model: Generic  512GB NVMe disk
* Primary Data Storage:  
    - Type: External USB HDD 
    - Capacity: 931.5 GB 
    - Model: 1TB EXTERNAL_USB (probablemente un Seagate o Western Digital) 
* RAM: 64 GB DDR5 (2 x 32 GB MRX5U560LKKD32G DIMMs)
    - Speed: 5600 MT/s (configurada a 5200 MT/s)
    - Type: Synchronous DDR5, Volatile memory

The operating system used was Ubuntu 24.04.2 LTS , the latest Long-Term Support (LTS) release at the time of analysis.
For most bioinformatic tools, we used Anaconda3 v2024-10-1 as the main environment manager.

## Step 1: Prepare the setup

Download and install Anaconda:
```bash
# Install basics

sudo apt update && sudo apt upgrade -y
sudo apt install build-essential
sudo apt install curl wget
sudo apt install git
sudo apt install vim nano
sudo apt install htop
sudo apt install zip unzip
sudo apt install default-jdk
sudo apt install perl

# Miniconda
cd ~
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh 
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
conda --version
python --version

# Chanels configuration
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
```
When prompted, select "Yes" to add Miniconda to the PATH automatically.


## Step 2: Install and activate bioinformatic tools
Once Minionda is installed, create the environment with the bioinformatic tools, Conda chanels priority, and dependencies for the raw reads analyses. 
Use the bioinfo-tools.yml file.

```bash
conda env create -f bioinfo-tools.yml

# Check the created environment
conda env list

# List the libraries in the
conda list environment
```
This will install only the packages mentioned.

## Step 3: Request data
The raw data can be downloaded from projects PRJNA1269778 and PRJNA954561.

```bash
# SRA Toolkit:
https://github.com/ncbi/sra-tools/wiki/Downloads
conda install -c bioconda sra-tools

# save the SRR.. files in sra_files.txt

# prefetch and fasterq-dump to download all the .sra files
cat sra_files.txt | while read id; do
    prefetch $id
    fasterq-dump $id --split-files --threads 10
done
```

Because sequencing documents are relatively large, pre-hash comparison is recommended to ensure that downloaded files are not corrupted.
