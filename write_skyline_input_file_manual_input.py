# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:53:01 2022

@author: nicola.freyer

Writes skyline inputn files from selected entries of subs file.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import re
import math
from datetime import datetime as dt

#%% Configuration

dir_path = "W:/Nicola/Data/fusA_mutants/220521_EFG_membrane_proteins/Band1_v220531/"

# Output file names:
filename_ssl = "skyline_input/skyline_input"
filename_fasta = "skyline_input/fasta_skyline"

mz_decimals = 1

# Create timestamp for output files:
now = dt.now()
timestamp = now.strftime("_v%Y-%m-%d_%H-%M-%S")

#%% Define functions

def fasta_line_short(row):
    entry = (">sp|addon|" + 
             str(row["protein"]) + "_" + 
             str(row["substitution"]).replace(" ", "_") + "_" + 
             str(row["Peptide_type"]) + "\n" +
             str(row["modified_sequence"]) + "\n")
    return entry

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

#%% Create target list from subs files

subs95 = pd.read_pickle(os.path.join(dir_path, "output/subs_unimod_95"))
subs80 = pd.read_pickle(os.path.join(dir_path, "output/subs_unimod_80"))

subs_diff = pd.merge(subs95, subs80, how="outer", indicator=True)
subs_diff = subs_diff[subs_diff["_merge"]=="right_only"]            # filters in entries from 80 that are not in 95

subs_diff = subs_diff[
    (subs_diff["mispairing"]!=False)                 # filters out non-cognates
    & (subs_diff["danger"]==False)                        # filters out all unimods
    # & (~subs_diff["substitution"].isin(subs_out))    # filters out selected unimods
    # & (subs_diff["protein"].str.match("tuf[AB]"))  # filters out non-EF-Tu peptides
    ]

# Clean data
subs_diff["mispairing"].fillna(0, inplace=True)
subs_diff["codon"].fillna("NNN", inplace=True)
subs_diff["position"].fillna(0, inplace=True)

#%% Write fasta file

# Remove replicate entries
columns_ex = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution", "Raw file"]
columns_gb = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution"]
subs_gb = subs_diff[columns_ex].groupby(columns_gb).count().reset_index()

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

pep_list = list(subs_diff["modified_sequence"].unique())

path_to_msmsscans = os.path.join(dir_path, "msmsScans.txt")
mms_iter = pd.read_csv(path_to_msmsscans, sep="\t", chunksize=10000, iterator=True, low_memory=False)
mms = pd.concat(chunk for chunk in mms_iter)
mms.reset_index(drop=True, inplace=True)

#%% Cross reference DP & BP
subs_diff["id"] = subs_diff.reset_index(drop=True).index
mms["id"] = mms.reset_index(drop=True).index

mms["m/z"] = mms["m/z"].round(mz_decimals)
subs_diff["m/z"] = subs_diff["m/z"].round(mz_decimals)

subs_merged = subs_diff.merge(mms, how="left", 
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

BP_list = list(subs_merged["DP Base Sequence"].unique())

mms_BPscans = mms[mms["Sequence"].isin(BP_list)]

output_BP = pd.DataFrame({
    "file": mms_BPscans["Raw file.raw"],
    "scan": mms_BPscans["Scan number"],
    "charge": mms_BPscans["Charge"],
    "sequence": mms_BPscans["Sequence"],
    "score": mms_BPscans["Score"],
    "modifications": mms_BPscans["Sequence"]
    })

output_BP = output_BP[output_BP["charge"] != 0]

# Diagnostics:
# output_BP.to_csv(os.path.join(dir_path, "Diagnostics/output_BP.csv"))
# gb1 = output_BP.groupby(["sequence"])["score"].count().reset_index()
# gb1.to_csv(os.path.join(dir_path, "Diagnostics/output_BP_gb1.csv"))

#%% Final output

output = pd.concat([output_BP, output_DP], ignore_index=True)

output["sequence"] = output["sequence"].str.replace("C", "C[+57.0]")
output["modifications"] = output["modifications"].str.replace("C", "C[+57.0]")

# Write ssl file
filename2 = dir_path + filename_ssl + timestamp + ".ssl"
output.to_csv(filename2, sep="\t", index=False)

# Diagnostic: 
output.to_csv(os.path.join(dir_path, "Diagnostics/output.csv"))
print("***INFO*** Number of rows in subs that have not been matched with mms row:")
print(output["score"].isna().sum())

