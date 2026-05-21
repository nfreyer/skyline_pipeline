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

#%% Variables

# Define logger colors
OKGREEN_TEXT = "\033[92m"
OKCYAN_TEXT = "\033[96m"
INFO_TEXT = "\033[93m"
ERROR_TEXT = "\033[91m"
ENDC_TEXT = "\033[0m"

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
        default = "subs_unimod.csv",
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

def log_me(text, log_file):
    """

    Parameters
    ----------
    text : str
        Text to be added to log document.

    Returns
    -------
    None.

    """
    timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")
    print(timestamp + "\t" + text, file=open(log_file, 'a'))

#%% Main

def main():
    time_start = dt.now()
    # Parser arguments
    args = get_arguments()
    dir_path = args.path
    subs_filename = args.filename
    mz_decimals = args.decimals
    
    # Check if directory exists
    if not os.path.isdir(dir_path):
        print(ERROR_TEXT + "ERROR ... The path specified does not exist." + ENDC_TEXT)
        sys.exit()
    
    # File names
    subs_path = dir_path + "output/" + subs_filename
    
    if not os.path.isfile(subs_path):
        print(ERROR_TEXT + "ERROR ... The file specified does not exist." + ENDC_TEXT)
        sys.exit()
    
    filename_fasta = "fasta_skyline"
    filename_ssl = "skyline_input"
    
    # Timestamp
    now = dt.now()
    timestamp = now.strftime("_v%Y-%m-%d_%H-%M-%S")
    timestamp1 = now.strftime("_v%Y%m%d")
    
    # Create output directory
    dir_path_out = os.path.join(dir_path, "skyline_input/")
    if os.path.isdir(dir_path_out):
        dir_input = input(OKCYAN_TEXT + "INPUT REQUIRED ... Output directory already exists. Enter y to create new: " + ENDC_TEXT)
        if dir_input.lower() == "y":
            print(INFO_TEXT + "INFO ... Overwriting output directory" + ENDC_TEXT)
            timestamp = dt.now().strftime("_v%Y-%m-%d_%H-%M-%S")
            dir_path_out = os.path.join(dir_path, "skyline_input" + timestamp + "/")
            os.mkdir(dir_path_out)
        else:
            print(ERROR_TEXT + "ERROR ... exit" + ENDC_TEXT)
            sys.exit()
    else:
        os.mkdir(dir_path_out)
        print(INFO_TEXT + "INFO ... Create output directory." + ENDC_TEXT)
    
    # Set up log file
    log_file = os.path.join(dir_path_out, "log.txt")
    with open(log_file, 'a') as file:
        file.write("CONFIGURATION\n")
        file.write("Input source:       MaxQuant\n")
        file.write(f"Directory:         {dir_path}\n")
        file.write(f"Input file name:   {subs_filename}\n")
        file.write(f"m/z accuracy:      {mz_decimals} decimal(s)\n")
    
    #%% Import subs, filter, clean
    log_me("### PROCESS STARTED ###", log_file = log_file)
    print(INFO_TEXT + "INFO ... Importing subs file." + ENDC_TEXT)
    
    subs = pd.read_csv(subs_path)
    
    # Dynamic subs filtering
    log_me("Started dynamic filtering of substitutions:", log_file = log_file)
    list_of_filters = []
    
    if args.mispairing:
        list_of_filters.append('mispairing != False')
        print(INFO_TEXT + "FILTERING ... Filter out non-cognate errors." + ENDC_TEXT)
        log_me("> Filter out non-cognate decoding errors", log_file = log_file)
    
    if args.danger:
        list_of_filters.append('danger == False')
        print(INFO_TEXT + "FILTERING ... Filter out unimods." + ENDC_TEXT)
        log_me("> Filter out unimods", log_file = log_file)
    
    if args.protein:
        protein_filter = 'protein.str.match("' + str(args.protein) + '")'
        list_of_filters.append(protein_filter)
        print(INFO_TEXT + 'FILTERING ... Filter in peptides from ' + str(args.protein) + '.' + ENDC_TEXT)
        log_me("> Filter in peptides from " + str(args.protein), log_file = log_file)
    
    if args.subs_in:
        subs_in = args.subs_in
        subs_in_listoflists = [list(x) for x in subs_in.split("+")]
        subs_in_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_in_listoflists]
        subs_in_filter = 'substitution in @subs_in_list'
        list_of_filters.append(subs_in_filter)
        print(INFO_TEXT + 'FILTERING ... Filter in substitutions: ' + ', '.join(subs_in_list) + ENDC_TEXT)
        log_me("> Filter in substitutions: " + ', '.join(subs_in_list), log_file = log_file)
        
    if args.subs_out:
        subs_out = args.subs_out
        subs_out_listoflists = [list(x) for x in subs_out.split("+")]
        subs_out_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_out_listoflists]
        subs_out_filter = 'substitution not in @subs_out_list'
        list_of_filters.append(subs_out_filter)
        print(INFO_TEXT +  'FILTERING ... Filter out substitutions: ' + ', '.join(subs_out_list) + ENDC_TEXT)
        log_me("> Filter out substitutions: " + ', '.join(subs_out_list), log_file = log_file)
    
    if args.free_text:
        ftf = args.free_text
        list_of_filters.append(ftf)
        print(INFO_TEXT + 'FILTERING ... Apply free text filter: ' + ftf + ENDC_TEXT)
        log_me("> Free text filter: " + ftf, log_file = log_file)
    
    query_cond = " & ".join(list_of_filters)
    
    if any([args.mispairing, args.danger, args.protein, args.subs_in, args.subs_out, args.free_text]):
        subs.query(query_cond, inplace=True, engine="python")
    
    if len(subs) == 0:
        print(ERROR_TEXT + 
              "ERROR ... No entries in subs after filtering, program was exited. Conflicting filters might have been applied." +
              ENDC_TEXT)
        log_me("No entries in subs after filtering. Exit.", log_file = log_file)
        sys.exit()    
    else:
        print(INFO_TEXT + "INFO ...", len(subs), "entries in subs after filtering." + ENDC_TEXT)
        log_me("> Finished filtering. " + str(len(subs)) + " entries left in subs", log_file = log_file)
        
        
    # Clean data
    subs["mispairing"] = subs["mispairing"].fillna(0)
    subs["codon"] = subs["codon"].fillna("NNN")
    subs["position"] = subs["position"].fillna(0)
       
    subs.to_csv(os.path.join(dir_path_out, "subs_filtered" + timestamp1 + ".csv"))
    # print("***INFO*** Nan Values in subs:")
    # print(subs.isna().sum())
    
    #%% Write fasta file
    
    print(INFO_TEXT + "INFO ... Preparing fasta file." + ENDC_TEXT)
    
    # Remove replicate entries
    columns_gb = ["protein", "DP Base Sequence", "modified_sequence", "position", "codon", "substitution"]
    subs_gb = subs[columns_gb + ["Raw file"]].groupby(columns_gb).count().reset_index()
    
    # Diagnostics: 
    # subs_gb.to_csv(os.path.join(dir_path_out, "subs_gb.csv"))
    
    subs_gb.rename({"modified_sequence": "Seq_DP", "DP Base Sequence": "Seq_BP"}, axis=1, inplace=True)
    subs_gb["id"] = subs_gb.index
    subs_long = pd.wide_to_long(subs_gb, stubnames="Seq_", i="id", j="Peptide_type", suffix=".+").reset_index(0, drop=True).reset_index(drop=False)
    subs_long.sort_values(by = ["protein", "position", "substitution", "Peptide_type"], inplace=True)
    subs_long.drop_duplicates(subset=["protein", "Peptide_type", "Seq_"], inplace=True, ignore_index=True)
    
    
    # subs_long["Num"] = subs_long.index
    subs_long["Num"] = np.nan
    subs_long.loc[subs_long["Peptide_type"] == "BP", "Num"] = subs_long[subs_long["Peptide_type"] == "BP"].reset_index().index + 1
    subs_long["Num"] = subs_long["Num"].ffill()
    
    subs_long["fasta"] = subs_long.apply(lambda row: fasta_line(row), axis=1)
    
    # Diagnostics:
    # subs_long.to_csv(os.path.join(dir_path_out, "subs_long.csv"))
    
    # Write fasta file
    filename1 = dir_path_out + filename_fasta + timestamp + ".fa"
    f = open(filename1, "a+")
    
    for index, row in subs_long.iterrows():
        f.write(row["fasta"])
    
    f.close()
    
    print(OKGREEN_TEXT + "INFO ... Exported fasta file successfully." + ENDC_TEXT)
    log_me("Exported fasta file with " + str(len(subs_long)) + " entries", log_file = log_file)
    
    #%% Import msms scans
    
    print(INFO_TEXT + "INFO ... Importing msms scans file." + ENDC_TEXT)
    
    path_to_msmsscans = os.path.join(dir_path, "msmsScans.txt")
    mms_iter = pd.read_csv(path_to_msmsscans, sep="\t", chunksize=10000, iterator=True, low_memory=False)
    mms = pd.concat(chunk for chunk in mms_iter)
    mms.reset_index(drop=True, inplace=True)
    
    log_me("Imported msms scans file.", log_file = log_file)
    
    #%% Cross reference DP & BP
    
    print(INFO_TEXT + "INFO ... Preparing ssl file." + ENDC_TEXT)
    
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
    # subs_merged.to_csv(os.path.join(dir_path_out, "subs_merged.csv"))
    # duplicates = subs_merged[subs_merged[["Raw file", "m/z", "Charge", "DP Base Sequence", "DP Probabilities"]].duplicated(keep=False)]
    # duplicates.to_csv(os.path.join(dir_path_out, "duplicates.csv"))
    
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
    # output_DP.to_csv(os.path.join(dir_path_out, "output_DP.csv"))
    
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
    # output_BP.to_csv(os.path.join(dir_path_out, "output_BP.csv"))
    # gb1 = output_BP.groupby(["sequence"])["score"].count().reset_index()
    # gb1.to_csv(os.path.join(dir_path_out, "output_BP_gb1.csv"))
    
    #%% Final output
    
    output = pd.concat([output_BP, output_DP], ignore_index=True)
    
    output["sequence"] = output["sequence"].str.replace("C", "C[+57.0]")
    output["modifications"] = output["modifications"].str.replace("C", "C[+57.0]")
    
    # Write ssl file
    print(INFO_TEXT + "INFO ... Writing ssl file." + ENDC_TEXT)
    
    filename2 = dir_path_out + filename_ssl + timestamp + ".ssl"
    output.to_csv(filename2, sep="\t", index=False)
    
    print(OKGREEN_TEXT + "INFO ... Exported ssl file successfully." + ENDC_TEXT)
    log_me("Exported ssl file", log_file = log_file)
    
    duration = dt.now() - time_start
    duration_str = str(duration).split(".")[0]
    
    # Diagnostics: 
    # output.to_csv(os.path.join(dir_path_out, "output.csv"))
    unmatched_rows = output["score"].isna().sum()
    if unmatched_rows == 0:
        print(OKGREEN_TEXT + "INFO ... Script finished. All subs rows have been matched with mms rows." + ENDC_TEXT)
        log_me(f"### PROCESS FINISHED - TOTAL RUNTIME: {duration_str} (HH:MM:SS) ###", log_file = log_file)
    else:
        print(ERROR_TEXT + "ERROR ... Number of rows in subs that have not been matched with mms row:", unmatched_rows, ENDC_TEXT)
        log_me(f"### ERROR ### {unmatched_rows} rows in subs have not been matched with mms rows", log_file = log_file)

#%%
    
if __name__ == '__main__':
    main()