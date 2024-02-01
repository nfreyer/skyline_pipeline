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
import sys
from datetime import datetime as dt

import argparse

os.system('color')

#%%

def get_arguments():
    """ Function to collect command line arguments.

    Returns
    -------
    args : argparse.Namespace
        The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description='Writes skyline input files from Pilpel script outputs.')

    parser.add_argument(
        "path",
        help = 'Directory path of input files.'
    )

    parser.add_argument(
        "-dec",
        "--decimals",
        type = int,
        default = 1,
        help = 'Number of decimals to use for m/z values when merging subs & msmsScans files. Default = 1'
    )
    
    parser.add_argument(
        "-fn",
        "--filename",
        default = "subs_unimod",
        help = 'Name of Pilpel script output file to use as subs input. Default = subs_unimod'
    )
    
    parser.add_argument(
        "-m",
        "--mispairing",
        action="store_true",
        help = 'Filters out non-cognate errors from subs file (mispairing =! False).'
    )
    
    parser.add_argument(
        "-d",
        "--danger",
        action="store_true",
        help = 'Filters out all unimods from subs file (danger == False).'
    )
    
    parser.add_argument(
        "-p",
        "--protein",
        help = 'Filters out all peptides not from the specified protein. Requires gene name / regex as input.'
    )

    parser.add_argument(
        "-si",
        "--subs_in",
        help = 'Filters IN all the specified amino acid substitutions (input format: "DE+ND+VIL" etc.).'
    )
    
    parser.add_argument(
        "-so",
        "--subs_out",
        help = 'Filters OUT all the specified amino acid substitutions (input format: "DE+ND+VIL" etc.).'
    )
    
    parser.add_argument(
        "-ft",
        "--free_text",
        help = 'Free text filter for additional filter condition(s). Multiple inputs need to be separated by "&". (input format: "column is in [list]", "column == regex".)'
    )

    return parser.parse_args()

def fasta_line(row):
    """
    Parameters
    ----------
    row : pandas DataFrame
        Single row of pandas DataFrame.

    Returns
    -------
    entry : String
        Entry for fasta file containing descriptor + sequence.

    """
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

def main():
    # Parser arguments
    args = get_arguments()
    dir_path = args.path
    
    if not os.path.isdir(dir_path):
        print("The path specified does not exist.")
        sys.exit()
    
    subs_filename = args.filename
    mz_decimals = args.decimals
    
    # File names
    subs_path = dir_path + "output/" + subs_filename
    
    if not os.path.isfile(subs_path):
        print("The file specified does not exist.")
        sys.exit()
    
    filename_fasta = "skyline_input/fasta_skyline"
    filename_ssl = "skyline_input/skyline_input"
    
    # Timestamp
    now = dt.now()
    timestamp = now.strftime("_v%Y-%m-%d_%H-%M-%S")
    
    #%%
    
    OKGREEN_TEXT = '\033[92m'
    INFO_TEXT = '\033[93m'
    ERROR_TEXT = '\033[91m'
    ENDC_TEXT = '\033[0m'
    
    #%% Import subs, filter, clean
    
    subs = pd.read_pickle(subs_path)
    
    # # Filters subs
    # subs_out = ["Q to E", "N to D", "F to Y", "T to E"] # unimods to exclude from analysis
    # subs = subs[
    #     (subs["mispairing"]!=False)                 # filters out non-cognates
    #     # & (subs["danger"]==False)                 # filters out all unimods
    #     # & (~subs["substitution"].isin(subs_out))  # filters out selected unimods
    #     # & (subs["protein"].str.match("tuf[AB]"))  # filters out non-EF-Tu peptides
    #     ]
    
    # Dynamic subs filtering
    list_of_filters = []
    
    if args.mispairing:
        list_of_filters.append('mispairing != False')
        print(INFO_TEXT + 
            "FILTERING ... Filter out non-cognate errors." + 
            ENDC_TEXT)
    
    if args.danger:
        list_of_filters.append('danger == False')
        print(INFO_TEXT +
            "FILTERING ... Filter out unimods." +
            ENDC_TEXT)
    
    if args.protein:
        protein_filter = 'protein.str.match("' + str(args.protein) + '")'
        list_of_filters.append(protein_filter)
        print(INFO_TEXT +
            'FILTERING ... Filter in peptides from ' + str(args.protein) + '.' +
            ENDC_TEXT)
    
    if args.subs_in:
        subs_in = args.subs_in
        subs_in_listoflists = [list(x) for x in subs_in.split("+")]
        subs_in_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_in_listoflists]
        subs_in_filter = 'substitution in @subs_in_list'
        list_of_filters.append(subs_in_filter)
        print(INFO_TEXT +
              'FILTERING ... Filter in substitutions: ' + ', '.join(subs_in_list) +
              ENDC_TEXT)
        
    if args.subs_out:
        subs_out = args.subs_out
        subs_out_listoflists = [list(x) for x in subs_out.split("+")]
        subs_out_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_out_listoflists]
        subs_out_filter = 'substitution not in @subs_out_list'
        list_of_filters.append(subs_out_filter)
        print(INFO_TEXT + 
              'FILTERING ... Filter out substitutions: ' + ', '.join(subs_out_list) +
              ENDC_TEXT)
    
    if args.free_text:
        ftf = args.free_text
        list_of_filters.append(ftf)
        print(INFO_TEXT + 
              'FILTERING ... Apply free text filter: ' + ftf + 
              ENDC_TEXT)
    
    query_cond = " & ".join(list_of_filters)
    
    if any([args.mispairing, args.danger, args.protein, args.subs_in, args.subs_out, args.free_text]):
        subs.query(query_cond, inplace=True, engine="python")
    
    if len(subs) == 0:
        print(ERROR_TEXT + 
              "ERROR ... No entries in subs after filtering, program was exited. Conflicting filters might have been applied." +
              ENDC_TEXT)
        sys.exit()    
    else:
        print(INFO_TEXT +
              "INFO ...", len(subs), "entries in subs after filtering." +
              ENDC_TEXT)
        
    # Clean data
    subs["mispairing"].fillna(0, inplace=True)
    subs["codon"].fillna("NNN", inplace=True)
    subs["position"].fillna(0, inplace=True)
    
    # Diagnostics:
    if os.path.isdir(os.path.join(dir_path, "Diagnostics/")):
        pass
    else:
        os.mkdir(os.path.join(dir_path, "Diagnostics/"))
        print(INFO_TEXT +
              "INFO ... Create diagnostics output directory." +
              ENDC_TEXT)
    
    subs.to_csv(os.path.join(dir_path, "Diagnostics/subs_filtered.csv"))
    # print("***INFO*** Nan Values in subs:")
    # print(subs.isna().sum())
    
    #%% Write fasta file
    
    # Remove replicate entries
    columns_gb = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution"]
    subs_gb = subs[columns_gb + ["Raw file"]].groupby(columns_gb).count().reset_index()
    
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
    # output.to_csv(os.path.join(dir_path, "Diagnostics/output.csv"))
    print(INFO_TEXT + 
          "INFO ... Number of rows in subs that have not been matched with mms row:", output["score"].isna().sum(),
          ENDC_TEXT)

if __name__ == '__main__':
    main()