#!/usr/bin/env python3
import os

# paths
aro_index = "aro_index.tsv"
fasta_files = "protein_fasta_protein_homolog_model.fasta"
output_dir = "~/databases/card/filt"

os.makedirs(output_dir, exist_ok = True)

# read index
index = {}

with open(aro_index) as f:
    for line in f:
        striped = line.strip()

        if not striped: # pass blank lines
            continue
        segments = striped.split("\t")

        if len(segments)<2: # pass when column < 2
            continue
        id_aro = segments[0] # link to fasta
        mechanism = segments[8] # column with ARO families
        index[id_aro] = mechanism
 
# read fasta
aro_mechs = set(index.values())

for mech_obj in aro_mechs:
    fasta_filt = []
    guardar = False
    seq_act = [] # blanck
    header = ""

    with open(fasta_files, "r") as f2:

        for line in f2:
            # save header if starts with >, is in index and is the same as mech_obj
            if line.startswith(">"):
                # Second loop, when guardar True and seq_act full, then save header and seq_act in fasta_filt
                if guardar and seq_act:
                    fasta_filt.append((header, "".join(seq_act)))

                striped2 = line.strip()
                segments2 = striped2.split("|")

                if len(segments2)>2:
                    aro_id = segments2[2]
                    
                    if aro_id in index and index[aro_id] == mech_obj:
                        guardar = True
                        header = striped2
                        seq_act = []
                    else:
                        guardar = False
                        seq_act = []
                else:
                    guardar = False
                    seq_act = []

            else:
                # save lines if them do not begin with > in seq_act
                if guardar:
                    seq_act.append(line.strip())

        # save header and seq_act when guardar True and seq_act full in the last loop
        if guardar and seq_act:
            fasta_filt.append((header, "".join(seq_act)))
    
    if fasta_filt:
        safe_name = mech_obj.replace(" ", "_").replace("/", "_")
        outname = os.path.join(output_dir, f"{safe_name}.fasta")

        with open(outname, "w")as fout:
            for header, seq_act in fasta_filt:
                fout.write(header + "\n")
                fout.write(seq_act+ "\n")





