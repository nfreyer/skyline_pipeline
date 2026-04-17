# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:53:01 2022

@author: nicola.freyer

Writes additional entries to add to automatically generated skyline input files.
"""

import pandas as pd
import os
from datetime import datetime as dt

#%% Configuration

dir_path = "W:/Nicola/Data/Error_cluster/220428_Double_errors/cleavage_pattern/"
protein = "tufA"

# Create timestamp for output files:
now = dt.now()
timestamp = now.strftime("_v%Y-%m-%d_%H-%M-%S")

#%% Define functions

def fasta_line_short(row):
    entry = (">sp|addon|" + 
             str(row["protein"]) + "_" + 
             str(row["substitution"]) + "_" + 
             str(row["Peptide_type"]) + "\n" +
             str(row["Peptide"]) + "\n")
    return entry

#%% Import peptide sequences

pep_df = pd.read_csv(os.path.join(dir_path, "peptide_list.csv"))
pep_df["protein"] = protein
pep_df["Peptide_type"] = "DP"
pep_df["substitution"] = pep_df["Protein Name"].apply(lambda x: x.split("|")[2])
pep_df["fasta"] = pep_df.apply(lambda row: fasta_line_short(row), axis=1)

#%% Write fasta file

filename = dir_path + "skyline_input/fasta_skyline_addon" + timestamp + ".fa"
f = open(filename, "a+")

for index, row in pep_df.iterrows():
    f.write(row["fasta"])

f.close()

#%% Import msms scans & filter for addon peptides

pep_list = list(pep_df["Peptide"].unique())

path_to_msmsscans = os.path.join(dir_path, "msmsScans.txt")
mms_iter = pd.read_csv(path_to_msmsscans, sep="\t", chunksize=10000, iterator=True)
mms = pd.concat(chunk for chunk in mms_iter)
mms.reset_index(drop=True, inplace=True)

mms_filtered = mms[mms["Sequence"].isin(pep_list)]

# Diagnostics:
mms_filtered.to_csv(os.path.join(dir_path, "Diagnostics/mms_filtered.csv"))

#%% Write output files

output = pd.DataFrame({
    "file": mms_filtered["Raw file"].apply(lambda x: x+".raw"),
    "scan": mms_filtered["Scan number"],
    "charge": mms_filtered["Charge"],
    "sequence": mms_filtered["Sequence"],
    "score": mms_filtered["Score"],
    "modifications": mms_filtered["Sequence"]
    })

print(output["charge"].eq(0).sum())
output = output[output["charge"] != 0]

output["sequence"] = output["sequence"].str.replace("C", "C[+57.0]")
output["modifications"] = output["modifications"].str.replace("C", "C[+57.0]")



# Diagnostics:
output.to_csv(os.path.join(dir_path, "Diagnostics/output.csv"))

# Write ssl file
output.to_csv(os.path.join(dir_path, "skyline_input/skyline_input_addon.ssl"), sep="\t", index=False, header=False)
