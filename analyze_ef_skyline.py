# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:05:30 2022

@author: nicola.freyer

Calculate the error frequencies from Skyline output data.
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

dir_path = "W:/Nicola/Data/Error_cluster/220428_Double_errors/"
filename = "error_frequency.csv"
filename_output = "Results/df_output.csv"
ident = "_quant"

#%% Import, filter & clean data

docgrid_head = pd.read_csv(os.path.join(dir_path, filename), nrows=10)
quants = [str(i) for i in docgrid_head.columns if ident in i]
docgrid_cols = quants + ["Protein Name", "Peptide"]
docgrid = pd.read_csv(os.path.join(dir_path, filename), usecols = docgrid_cols)

# Diagnostic
print("Number of nan values in input data: " + str(docgrid.isna().sum().sum()))

# Clean data
docgrid.fillna(0, inplace=True)

# Rename columns
quants_new = ["stub " + x.split(" ")[0] for x in quants]
quants_dict = {quants[i]: quants_new[i] for i in range(len(quants))}
docgrid.rename(quants_dict, axis=1, inplace=True, errors="raise")

# New columns
docgrid["id"] = docgrid.index
docgrid["pep_id"] = docgrid["Protein Name"].apply(lambda x: x.split("|")[1])
docgrid["type"] = docgrid["Protein Name"].apply(lambda x: x.split("|")[-1]).apply(lambda y: y.split("_")[-1])


#%% Reformat data & calculate EF
                                                                                                                                                                                                                                                                                  
docgrid_long = pd.wide_to_long(docgrid, stubnames = "stub ", 
                               i = "id", j = "Sample", suffix=".+").reset_index(1).reset_index(drop=True)
docgrid_long.rename({"stub ": "Peak area"}, axis=1, inplace=True)

# Diagnostic
# docgrid_long.to_csv(os.path.join(dir_path, "Results/docgrid_long.csv"))

df_BP = docgrid_long.loc[docgrid_long["type"]=="BP"]
df_DP = docgrid_long.loc[docgrid_long["type"]=="DP"]

df_merged = df_DP.merge(df_BP, how="left", on = ["Sample", "pep_id"], suffixes = ("_DP", "_BP"))


df_merged["EF"] = df_merged["Peak area_DP"] / (df_merged["Peak area_DP"] + df_merged["Peak area_BP"])
df_merged["Peptide_BP"].fillna("nan", inplace=True)

# Diagnostic
# df_merged.to_csv(os.path.join(dir_path, "Results/df_merged.csv"))
print("Number of DPs lost during merge: " + str(len(df_DP)-len(df_merged)))

#%% Pivot table to tidy table format & imputate

index_list = ["pep_id", "Protein Name_DP", "Peptide_BP", "Peptide_DP"]
df_tt = pd.pivot_table(df_merged, values = 'EF', index = index_list, 
                       columns = 'Sample', aggfunc = "mean").reset_index()

samples = [str(i) for i in df_tt.columns if ident in i]
df_tt[samples] = df_tt[samples].where(df_tt[samples] >= 0.000001, 0.000001, errors="raise")

# Diagnostic
# df_tt.to_csv(os.path.join(dir_path, "Results/df_tt.csv"))
print("Number of DPs lost during pivot: " + str(len(df_DP["Protein Name"].unique())-len(df_tt)))

#%% Filter for regulation

gb_list = [samples[i:i+3] for i in range(0, len(samples), 3)]
samples_short = ["_".join(x.split("_")[0:-1]) for x in samples]
# samples_unique = list(set(samples_short))
samples_unique = list(dict.fromkeys(samples_short))
samples_dict = {samples_unique[i]: gb_list[i] for i in range(len(samples_unique))}
# samples_tup = [(samples_short[i], samples[i]) for i in range(len(samples))]

df_median = df_tt.loc[:,["pep_id", "Protein Name_DP", "Peptide_BP", "Peptide_DP"]]
df_median[samples_unique] = np.nan

for key, value in samples_dict.items():
    df_median[key] = df_tt.loc[:,value].median(axis=1)

df_median["Reg_Mut"] = df_median["eng_P610L_256Apr"] / df_median["P610L_0Apr"]
df_median["Reg_WT"] = df_median["engwt_16"] / df_median["eng_wt_0Apr"]

# Diagnostic
# df_median.to_csv(os.path.join(dir_path, "Results/df_median.csv"))
# print(df_median)

df_filtered = df_median[df_median["Reg_Mut"] > 2].copy(deep=True).reset_index(drop=True)

# Diagnostic
print("Number of DP before regulation filter: " + str(len(df_median)))
print("Number of DP after regulation filter: " + str(len(df_filtered)))
# df_filtered.to_csv(os.path.join(dir_path, "Results/df_filtered.csv"))

print()

#%% Add addon EFs, pivot, filter for regulation

# Calculate median BP peak area per sample
id_list = list(df_filtered["pep_id"].unique())
medianEF_series = df_BP[df_BP["pep_id"].isin(id_list)].groupby("Sample")["Peak area"].median()
medianEF_dict = medianEF_series.to_dict()

# Fill BP peak area for addon peptides with median of all BPs per sample
df_merged.loc[df_merged["pep_id"]=="addon", "Peak area_BP"] = df_merged.loc[df_merged["pep_id"]=="addon", "Sample"].map(medianEF_dict)
df_addon = df_merged.loc[df_merged["pep_id"]=="addon", :].copy(deep=True)

# Calculate EF
df_addon["EF"] = df_addon["Peak area_DP"] / (df_addon["Peak area_DP"] + df_addon["Peak area_BP"])

# Pivot table to tidy table format
df_addon_tt = pd.pivot_table(df_addon, values = 'EF', index = index_list, 
                             columns = 'Sample', aggfunc = "mean").reset_index()

# Imputation
df_addon_tt[samples] = df_addon_tt[samples].where(df_addon_tt[samples] >= 0.000001, 0.000001, errors="raise")

# Calculate median
df_addon_median = df_addon_tt.loc[:,["pep_id", "Protein Name_DP", "Peptide_BP", "Peptide_DP"]]
df_addon_median[samples_unique] = np.nan

for key, value in samples_dict.items():
    df_addon_median[key] = df_addon_tt.loc[:,value].median(axis=1)

df_addon_median["Reg_Mut"] = df_addon_median["eng_P610L_256Apr"] / df_addon_median["P610L_0Apr"]
df_addon_median["Reg_WT"] = df_addon_median["engwt_16"] / df_addon_median["eng_wt_0Apr"]

df_addon_filtered = df_addon_median[df_addon_median["Reg_Mut"] > 2].copy(deep=True).reset_index(drop=True)

# Diagnostic
print("Number of BPs used for median: " + str(len(id_list)))
print("Number of addon DPs: " + str(len(df_addon_tt)))
print("Number of regulated addon DPs: " + str(len(df_addon_filtered)))

# %% Combine results files and export

df_output = pd.concat([df_filtered, df_addon_filtered], ignore_index=True)
# df_output = df_filtered

df_output["protein"] = df_output["Protein Name_DP"].apply(lambda x: x.split("|")[-1]).apply(lambda y: y.split("_")[0])

# Diagnostic
df_output.to_csv(os.path.join(dir_path, filename_output))
# print(df_filtered._is_copy)

#%% Selective deletion

# df_output = df_output[df_output["Protein Name_DP"] != "sp|30|tufA_237_K_to_N_AAA_DP"]

# Diagnostic
# df_output.to_csv(os.path.join(dir_path, filename_output))

#%% Plot results

colors = ["#C70039", "#9C9C9C", "#22B241", "#3A78D8"]
p = sns.color_palette(colors)

f, ax = plt.subplots(figsize=(5,5))
ax.axline((1, 1), slope=1, color="k", linestyle="--")

g = sns.scatterplot(data=df_output, x="engwt_16", y="eng_P610L_256Apr", hue="protein", 
                    # marker="o", color ="grey", 
                    edgecolor="k", palette=p)

g.set(xscale="log")
g.set(yscale="log")
plt.ylim(0.000001, 1)
plt.xlim(0.000001, 1)
plt.show()
f.savefig(os.path.join(dir_path, "Results/data.png"), dpi=200)

