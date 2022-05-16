# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:14:35 2022

@author: nicola.freyer

Writes skyline input files (.ssl + .fasta) from outputs of MaxQuant search (msmsScans) 
and Pilpel script (subs_unimods).
"""

import pandas as pd
import numpy as np
import os
from datetime import datetime as dt

#%% Configuration

# dir_path = "W:/Nicola/Data/Method_optimization/220129_EFG_HighEF_ident_10ppm/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/EF-Tu/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/protein1/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/protein2/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/protein3/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/protein4/"
# dir_path = "W:/Nicola/Data/Error_cluster/220427_Re-analysis/proteinH/"
dir_path = "W:/Nicola/Data/Error_cluster/220428_Double_errors/dependent_peptides/"

# Input files
subs_filename = "subs_unimod"
subs_path = dir_path + "output/" + subs_filename

mz_decimals = 1

# Output files
filename_fasta = "skyline_input/fasta_skyline_unimods_all"
filename_ssl = "skyline_input/skyline_input_unimods_all"

# Timestamp
now = dt.now()
timestamp = now.strftime("_v%Y-%m-%d_%H-%M-%S")

#%% Functions

def fasta_line(row):
    if row["Peptide_type"] == "DP":
        entry = (">sp|" + str(int(row["Num"])) + "|" + 
                 str(row["protein"]) + "_" + 
                 str(int(row["position"])) + "_" + 
                 row["substitution"].replace(" ", "_")  + "_" + 
                 str(row["codon"]) + "_" + 
                 str(row["Peptide_type"]) + "\n" + 
                 str(row["Seq_"]) + "\n")
    else:
        entry = (">sp|" + str(int(row["Num"])) + "|" + 
                 str(row["protein"]) + "_" + 
                 str(int(row["position"])) + "_" + 
                 str(row["Peptide_type"]) + "\n" + 
                 str(row["Seq_"]) + "\n")
    return entry

#%% Import subs, filter, clean

subs = pd.read_pickle(subs_path)

# Filters subs
subs_out = ["Q to E", "N to D", "F to Y", "T to E"] # unimods to exclude from analysis
subs = subs[
    (subs["mispairing"]!=False)                 # filters out non-cognates
    # & (subs["danger"]==False)                 # filters out all unimods
    & (~subs["substitution"].isin(subs_out))    # filters out selected unimods
    # & (subs["protein"].str.match("tuf[AB]"))  # filters out non-EF-Tu peptides
    ]

# Clean data
subs["mispairing"].fillna(0, inplace=True)
subs["codon"].fillna("NNN", inplace=True)
subs["position"].fillna(0, inplace=True)

# Diagnostics:
# subs.to_csv(os.path.join(dir_path, "Diagnostics/subs_cleaned.csv"))
print("***INFO*** Nan Values in subs:")
print(subs.isna().sum())

#%% Write fasta file

# Remove replicate entries
columns_ex = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution", "Raw file"]
columns_gb = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution"]
subs_gb = subs[columns_ex].groupby(columns_gb).count().reset_index()

# Diagnostics: 
# subs_gb.to_csv(os.path.join(dir_path, "Diagnostics/subs_gb.csv"))

subs_gb.rename({"modified_sequence": "Seq_DP", "DP Base Sequence": "Seq_BP"}, axis=1, inplace=True)
subs_gb["id"] = subs_gb.index
subs_long = pd.wide_to_long(subs_gb, stubnames="Seq_", i="id", j="Peptide_type", suffix=".+").reset_index(0, drop=True).reset_index(drop=False)
subs_long.sort_values(by = ["protein", "position", "substitution", "Peptide_type"], inplace=True)
subs_long.drop_duplicates(subset=["protein", "Peptide_type", "Seq_"], inplace=True, ignore_index=True)


# subs_long["Num"] = subs_long.index
subs_long["Num"] = np.nan
subs_long.loc[subs_long["Peptide_type"] == "BP", "Num"] = subs_long[subs_long["Peptide_type"] == "BP"].reset_index().index + 1
subs_long["Num"].fillna(method="ffill", inplace=True)

subs_long["fasta"] = subs_long.apply(lambda row: fasta_line(row), axis=1)

# Diagnostics:
# subs_long.to_csv(os.path.join(dir_path, "Diagnostics/subs_long.csv"))

# Write fasta file
filename1 = dir_path + filename_fasta + timestamp + ".fa"
f = open(filename1, "a+")

for index, row in subs_long.iterrows():
    f.write(row["fasta"])

f.close()

#%% Import msms scans

path_to_msmsscans = os.path.join(dir_path, "msmsScans.txt")
mms_iter = pd.read_csv(path_to_msmsscans, sep="\t", chunksize=10000, iterator=True, low_memory=False)
mms = pd.concat(chunk for chunk in mms_iter)
mms.reset_index(drop=True, inplace=True)

#%% Cross reference DP & BP
subs["id"] = subs.reset_index(drop=True).index
mms["id"] = mms.reset_index(drop=True).index

mms["m/z"] = mms["m/z"].round(mz_decimals)
subs["m/z"] = subs["m/z"].round(mz_decimals)

subs_merged = subs.merge(mms, how="left", 
           left_on = ["Raw file", "m/z", "Charge", "DP Base Sequence", "DP Probabilities"], 
           right_on = ["Raw file", "m/z", "Charge", "DP base sequence", "DP probabilities"],
           validate = "many_to_many"
           )

# Diagnostics:
# subs_merged.to_csv(os.path.join(dir_path, "Diagnostics/subs_merged.csv"))
# duplicates = subs_merged[subs_merged[["Raw file", "m/z", "Charge", "DP Base Sequence", "DP Probabilities"]].duplicated(keep=False)]
# duplicates.to_csv(os.path.join(dir_path, "Diagnostics/duplicates.csv"))


#%% DP output table

subs_merged["Raw file.raw"] = subs_merged["Raw file"] + ".raw"

output_DP = pd.DataFrame({
    "file": subs_merged["Raw file.raw"],
    "scan": subs_merged["Scan number"],
    "charge": subs_merged["Charge"],
    "sequence": subs_merged["modified_sequence"],
    "score": subs_merged["DP score"],
    "modifications": subs_merged["modified_sequence"]
    })

# Diagnostics:
# output_DP.to_csv(os.path.join(dir_path, "Diagnostics/output_DP.csv"))

#%% BP output table

mms["Raw file.raw"] = mms["Raw file"] + ".raw"

gb = subs_merged.groupby(["DP base scan number", "Raw file"])["m/z"].count().reset_index()

# Create identifier tuple consisting of scan number & raw file name (without .raw)
ident_list = list(zip(gb["DP base scan number"], gb["Raw file"]))
mms["ident_tuple"] = list(zip(mms["Scan number"], mms["Raw file"]))

mms_BPscans = mms[mms["ident_tuple"].isin(ident_list)]

# Diagnostics:
# mms_BPscans.to_csv(os.path.join(dir_path, "Diagnostics/mms_BPscans.csv"))

output_BP = pd.DataFrame({
    "file": mms_BPscans["Raw file.raw"],
    "scan": mms_BPscans["Scan number"],
    "charge": mms_BPscans["Charge"],
    "sequence": mms_BPscans["Sequence"],
    "score": mms_BPscans["Score"],
    "modifications": mms_BPscans["Sequence"]
    })

# Diagnostics:
# output_BP.to_csv(os.path.join(dir_path, "Diagnostics/output_BP.csv"))

#%% Final output

output = pd.concat([output_BP, output_DP], ignore_index=True)

output["sequence"] = output["sequence"].str.replace("C", "C[+57.0]")
output["modifications"] = output["modifications"].str.replace("C", "C[+57.0]")

# Write ssl file
filename2 = dir_path + filename_ssl + timestamp + ".ssl"
output.to_csv(filename2, sep="\t", index=False)

# Diagnostic: 
# output.to_csv(os.path.join(dir_path, "Diagnostics/output.csv"))
print("***INFO*** Number of rows in subs that have not been matched with mms row:")
print(output["score"].isna().sum())
