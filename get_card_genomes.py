#!/usr/bin/env python3
import csv
import time
import re
from Bio import Entrez
from Bio import SeqIO
from io import StringIO
import os


# Config
Entrez.email = "direction@mail.com"
NCBI_API_KEY = "..."
BATCH_SIZE = 100 
OUTPUT_FILE = "~/databases/card/taxonomic_info.tsv"


def extract_protein_accessions(tsv_file):
    accessions = set()
    with open(tsv_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        
        for row in reader:
            if len(row) > 6:
                acc_field = row[6].strip()
                if acc_field and acc_field != "NA":
                    for acc in re.split(r'[,\s]+', acc_field):
                        if acc.startswith(("NP_", "XP_", "WP_", "YP_")):
                            accessions.add(acc)
    return accessions


def get_taxonomy_info(accessions):
    with open(OUTPUT_FILE, "w") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["Protein_Accession", "Organism", "TaxID", "Source_Accession"])
        
        for i in range(0, len(accessions), BATCH_SIZE):
            batch = list(accessions)[i:i+BATCH_SIZE]
            attempt = 0
            
            while attempt < 3:
                try:
                    handle = Entrez.efetch(
                        db="protein",
                        id=",".join(batch),
                        rettype="gb",
                        retmode="text",
                        api_key=NCBI_API_KEY
                    )
                    records = SeqIO.parse(StringIO(handle.read()), "genbank")
                    
                    for record in records:
                        acc = record.id.split(".")[0]
                        if "organism" in record.annotations:
                            org = ' '.join(record.annotations["organism"].split()[:2])
                            taxid = record.annotations.get("taxid", [""])[0]
                            source_acc = record.id
                            writer.writerow([acc, org, taxid, source_acc])
                    
                    break
                    
                except Exception as e:
                    print(f"Error in batch {i//BATCH_SIZE + 1}: {str(e)}")
                    attempt += 1
                    time.sleep(15)
            
            time.sleep(1)


def main():
    # Extract CARD accesion
    accessions = extract_protein_accessions("aro_index.tsv")
    print(f"Found {len(accessions)}")
    
    # Save taxonomic info
    get_taxonomy_info(accessions)
    print(f"\nSaved in: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()