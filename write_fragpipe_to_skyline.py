# -*- coding: utf-8 -*-
"""
Created on Thu Feb 12 13:54:12 2026

formats fragpipe (msfragger, mass-offset search) results to subs output format
to allow feeding into skyline_input script

Requires config file
Requires subs file for each codon table

@author: nicola.freyer
"""

import pandas as pd
import numpy as np
import os
import sys
import re

from Bio import SeqIO # only needed if codon information is required #TODO
import Bio.Data.CodonTable as cd
from Bio.Seq import Seq

from functools import partial
from datetime import datetime as dt
import argparse
import yaml

os.system('color')

positional_probability_cutoff = 0.95 # minimal threshold for delta mass localization probability # TODO
# literally where is this in fragpipe

#%% Variables

# Define logger colors
OKGREEN_TEXT = "\033[92m"
OKCYAN_TEXT = "\033[96m"
INFO_TEXT = "\033[93m"
ERROR_TEXT = "\033[91m"
ENDC_TEXT = "\033[0m"

#%% Functions

def hamming(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def get_codon_table(codons, amino_acids):
    return dict(list(zip(codons, amino_acids)))

def get_inverted_codon_table(ct):
    inv_codon_table = {}
    for k, v in ct.items():
        inv_codon_table[v] = inv_codon_table.get(v, [])
        inv_codon_table[v].append(k)
    return inv_codon_table

def codonify(seq):
    """
    Turns nucleotide sequence into list of triplet codons

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    list
        list of triplet codons.

    """
    seq = str(seq)
    l = len(seq)
    return [seq[i : i + 3] for i in range(0, l, 3)]

def is_subs(row, tol, df):
    """
    Takes a row of the psm.tsv style file (df) as input and checks this row's 
    corresponding origin aa ("aa_origin" column) and modification delta mass
    ("deltaM" column) against the reference table for aa substitutions 
    to identify delta mass & aa combinations that could be produced by 
    an amino acid substitution.

    Parameters
    ----------
    row : pandas row
        row of the psm.tsv file that is currently tested.
    df : pandas table
        reference table of aa subsitutions and corresponding delta masses
        contains "delta_mass" and "aa_origin" columns
    
    Returns
    -------
    result : str (or list of str)
        empty if no aa subs is possible, 
        list of subs types if multiple subs are possible (impossible with mass tol of 0.005)
        str for single possible subs type

    """
    aa_origin = row["aa_origin"]
    delta_mass = row["deltaM"]
    dmass_min = delta_mass - tol
    dmass_max = delta_mass + tol
    hit = df.loc[ (df["delta_mass"].between(dmass_min, dmass_max)) & 
                     (df["aa_origin"] == aa_origin) ]
    if hit.empty:
        result = ""
    elif len(hit) > 1:
        # result = hit.loc[:,"Sub"].to_list() # should not be possible, remove? will cause issues downstream
        result = "ERROR, multiple hits"
    else:
        result = hit["Sub"].to_list()[0]
    
    return result

def is_gene(record, ct_id):
    """
    Checks if a record entry is a gene (starts with start codon, has triplet codons, and ends with stop codon)

    Parameters
    ----------
    record : SeqIO record
        fasta file record.
    ct_id : int
        NCBI codon table.

    Returns
    -------
    bool
        True if the record fullfills all three criteria to qualify as a gene, False otherwise.

    """
    if len(record.seq) % 3 != 0:
        return False
    if not record.seq[:3] in set(cd.generic_by_id[int(ct_id)].start_codons):
        return False
    if record.seq[-3:].translate() != "*":
        return False
    return True

def is_mispairing(row, mask=None, mis_dict=None):
    """
    Returns whether the substitution is mispairing or misloading, based on the
    near-cognate mask.
    """
    codon = row["codon"]
    destination = row["aa_dest"]
    if pd.notnull(codon) and pd.notnull(destination):
        if (codon in mask.index) and destination:
            return mask.loc[codon, destination] # TODO # throws problem
        else: # if codon is not in the mask index # does not control for cases where destination is empy, but when would it even be empty?
            return mis_dict[row["is_sub"]] # map mispairing dict from subs
    else:
        return float("NaN")

def get_arguments():
    """ 
    Function to collect command line arguments.
    Adapted from write_skyline_input_file.py

    Returns
    -------
    args : argparse.Namespace
        The parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description='Writes skyline input files from fragpipe (msfragger, mass-offset search) output.')

    parser.add_argument(
        "path",
        help = 'Directory path of input files.'
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

def read_yaml(config_path):
    """
    Reads in yaml file and outputs a dictionary with its contents.
    
    Parameters
    ----------
    config_path : str
        path to config.yaml file including the file name.

    Returns
    -------
    dict
        Dictionary containing the input of the yaml config file.
    """
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def fasta_line(row):
    """
    Function to write a fasta file style entry from a DataFrame style table.
    Adapted from write_skyline_input_file.py (adjusted column names)
    
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
                 str(row["Gene"]) + "_" +                       # protein name
                 str(int(row["pos"])) + "_" +                   # subs pos in protein
                 row["Sub"].replace(" ", "_")  + "_" +          # subs type (no whitespace)
                 str(row["Peptide_type"]) + "\n" +              # removed codon information
                 str(row["Seq_"]) + "\n")
    else:
        entry = (">sp|" + str(int(row["Num"])) + "|" + 
                 str(row["Gene"]) + "_" +                       # protein name
                 str(int(row["pos"])) + "_" +                   # subs pos in protein
                 str(row["Peptide_type"]) + "\n" + 
                 str(row["Seq_"]) + "\n")
    return entry

#%% Main

def main():
    # Parse arguments
    args = get_arguments()
    dir_path = args.path

    # Check if directory exists
    if not os.path.isdir(dir_path):
        print("The path specified does not exist.")
        sys.exit()

    # Checks if the configuration file exists.
    config_path = os.path.join(dir_path, "config.yaml")
    
    if not os.path.isfile(config_path):
        print(ERROR_TEXT + "The configuration file specified does not exist." + ENDC_TEXT)
        sys.exit()
    
    # Import values from configuration file
    config = read_yaml(config_path)
    tol = float(config["tol"])
    ct_id = int(config["codon_table"])
    cds_path = config["path_to_cds"]

    # Create output directory
    dir_path_out = os.path.join(dir_path, "output_fragpipe/")
    if os.path.isdir(dir_path_out):
        dir_input = input(OKCYAN_TEXT + "INPUT REQUIRED ... Output directory already exists. Press y to create new: " + ENDC_TEXT)
        if dir_input.lower() == "y":
            print(INFO_TEXT + "INFO ... Overwriting output directory" + ENDC_TEXT)
            timestamp = dt.now().strftime("_v%Y-%m-%d_%H-%M-%S")
            dir_path_out = os.path.join(dir_path, "output_fragpipe" + timestamp + "/")
            os.mkdir(dir_path_out)
        else:
            print(ERROR_TEXT + "ERROR ... exit" + ENDC_TEXT)
            sys.exit()
    else:
        os.mkdir(os.path.join(dir_path, "output_fragpipe/"))
        print(INFO_TEXT +
              "INFO ... Create output directory." +
              ENDC_TEXT)
            
    #%% Import files & create psm sum file
    print(INFO_TEXT + 
          "INFO ... Importing fragpipe results.",
          ENDC_TEXT)
    
    with open(os.path.join(dir_path, "filelist_proteinprophet.txt"), "r") as file:
        filenames = [line.rstrip().split("\\")[-2] for line in file]
    
    print(INFO_TEXT + 
          "INFO ... Raw files considered for analysis: \n \t - " +
          "\n \t - ".join(filenames) + 
          ENDC_TEXT)
    
    # reads in all psm.tsv files (generated seperately for each raw file) & concats into one df
    df = pd.DataFrame()
    for file in filenames:
        dfiter = pd.read_csv(os.path.join(dir_path,file,'psm.tsv'), sep='\t', low_memory=False)
        df = pd.concat([df,dfiter])
    df.reset_index(drop=True,inplace=True)

    # Prepare df for ssl file
    df["Scan"] = [x.split(".")[1] for x in df["Spectrum"]]
    df["File_raw"] = [x.split(".")[0] + ".raw" for x in df["Spectrum"]]
    df["RT_min"] = [float(x) / 60 for x in df["Retention"]]
    
    df["mods_list"] = [x.split(", ") if type(x) == str else "" for x in df["Assigned Modifications"] ]
    df["mods_clean"] = df["mods_list"].map(lambda x: [a for a in x if r"C(57.021" not in a])       # filters out CamCys from mods
    df["mods_clean2"] = df["mods_clean"].map(lambda x: [a for a in x if r"N-term(" not in a])       # filters out Nterm mods from mods

    # Diagnostic
    df.to_csv(os.path.join(dir_path_out, "psm_sum" + dt.now().strftime("_v%Y%m%d") + ".csv"))    
    
    print(OKGREEN_TEXT + 
          "INFO ... Exported psm_sum file successfully.",
          ENDC_TEXT)
    
    #%% Read in external input
    
    # TODO make modular: generate subs table for each codon table externally and read in? 
    # or auto-generate for each and give option to read in manually curated one?
    
    # Substituion matrix, generated by "generate_mastertable.py" plus manually curated for additional artifacts
    # Pandas format; columns for Sub ("X to Y"), aa_origin ("X"), aa_dest ("Y"), delta_mass (float), mispairing (bool), danger (bool)
    subs_table = pd.read_csv("W:/Nicola/Scripts/my_scripts/skyline_pipeline/Substitution_matrix_v241104.csv", index_col = 0)

    dict_mispairing = {k:v for k, v in zip(subs_table["Sub"],subs_table["mispairing"])}
    dict_danger     = {k:v for k, v in zip(subs_table["Sub"],subs_table["danger"])}
    
    #%% Import cDNA file
    
    print(INFO_TEXT + 
          "INFO ... Importing cDNA file.",
          ENDC_TEXT)
    
    # From Mordret substitution script
    record_dict = {}

    for record in SeqIO.parse(open(cds_path, "r"), "fasta"):
        record.seq = record.seq.upper()
        if is_gene(record, ct_id = ct_id):
            bits = record.description.split(" ")
            for i in bits:
                if "gene_symbol" in i:
                    record.name = i.split(":")[-1]
            record_dict[record.name] = record

    #%% filters df for psm where additional delta mass was detected
    
    print(INFO_TEXT + 
          "INFO ... Identify amino acid substitutions.",
          ENDC_TEXT)
    
    # Set up variables for misreading assessment
    bases = "TCAG"
    codons = [a + b + c for a in bases for b in bases for c in bases]
    aas = [str(Seq(x).translate(table=ct_id)) for x in codons]
    amino_acids = "".join(aas)
    codon_table = get_codon_table(codons, amino_acids)
    inverted_codon_table = get_inverted_codon_table(codon_table)
    inverted_codon_table["L"] = inverted_codon_table["L"] + inverted_codon_table["I"]

    # Create a mask to assess mispairing based on the nucleotide codon & the destination amino acid
    mask = pd.DataFrame(data=False, index=codons, columns=list("ACDEFGHKLMNPQRSTVWY"), dtype=bool)

    for label in codons:
        # True/False: Does this codon have exactly one base pair mismatch to all possible codons (in order)?
        near_cognates = np.array([hamming(i, label) == 1 for i in codons]) 
        # List unique amino acids corresponding to near cognate codons (check in order, only report set)
        reachable_aa = set(np.array(list(amino_acids))[near_cognates]) 
        # True/False: Is the amino acid reachable from the original codon through a single base pair mismatch?
        mask.loc[label] = [i in reachable_aa for i in "ACDEFGHKLMNPQRSTVWY"] 
    
    # Removes "near cognates" that produce the same amino acid (silent)
    # Seems to be redundant from testing? not sure, leave in to be on the safe side for other codon tables
    for label in mask.index:
        for col in mask.columns:
            if label in inverted_codon_table[col]:
                mask.loc[label, col] = False

    # Create new df with specifically psms containing a delta mass
    # copy() to avoid SettingWithCopyWarning 
    # dropna redudndant??
    df_subs = df.copy().dropna(subset=["Assigned Modifications"]) 
    
    ### TODO: A LOT OF POTENTIAL HERE: 
    
    df_subs = df_subs.loc[[len(x) == 1 for x in df_subs["mods_clean2"]]] # removes double modifications
    
    ###
    
    # find position of modificiation using "\d" (digit) and positive lookahead for "\w" (single letter) and "(" open bracket
    df_subs["pos_peptide"] = [re.search(r"\d+(?=\w\()", x[0]).group(0) for x in df_subs["mods_clean2"]]
    # add position in peptide and peptide start for mod position in protein
    # currently using amino acid numbering (1st AA = 1), not pythonic indexing (1st AA = 0)
    df_subs["pos"] = [int(x)+int(y)-1 for x, y in zip(df_subs["pos_peptide"], df_subs["Protein Start"])]
    # finds modified aa in mod syntax
    df_subs["aa_origin"] = [re.search(r"(?<=\d)[A-Z](?=\()", x[0]).group(0) for x in df_subs["mods_clean2"]]
    # finds delta mass in mod syntax
    df_subs["deltaM"] = [re.search(r"(?<=\()\S+(?=\))", x[0]).group(0) for x in df_subs["mods_clean2"]]
    # transforms delta mass from str to float
    df_subs["deltaM"] = df_subs["deltaM"].astype(float)
    
    # Identifies triplet codon from cDNA file
    df_subs["codon"] = [str(record_dict[g].seq[(p*3)-3:(p*3)])
                        if g in record_dict.keys() else "NNN"   # for unassigned genes (e.g. contaminations)
                        for g,p in zip(df_subs["Gene"], df_subs["pos"])]
    
    # partial function to apply is_subs function with default ref table
    is_subs_p = partial(is_subs, tol = tol, df=subs_table)
    # identify potential aa subs from fragpipe-identified modfication
    df_subs["is_sub"] = df_subs.apply(is_subs_p, 1)
    
    # Creates new df with subs only
    df_DP = df_subs[df_subs["is_sub"]!=""].copy()
    
    # Handles Ile/Leu ambiguity
    df_DP["Sub"] = df_DP["is_sub"].replace({"to L": "to I/L"}, regex=True)
    # finds destination aa from identified sub
    df_DP["aa_dest"] = [x.split(" ")[-1] for x in df_DP["is_sub"]]
    # maps danger boolean from ref table
    df_DP["danger"] = df_DP["is_sub"].map(dict_danger, na_action="ignore")
    
    # Assess mispairing
    
    # old version (considers only amino acids):
    # df_subs["mispairing"] = df_subs["is_sub"].map(dict_mispairing, na_action="ignore")
    
    # new version (considers codon identity too):
    is_mispairing_p = partial(is_mispairing, mask=mask, mis_dict=dict_mispairing)
    df_DP["mispairing"] = df_DP.apply(is_mispairing_p, 1)   
    
    # Writes modified sequence of DP
    df_DP["modified_sequence"] = [x[:int(y)-1] + z + x[int(y):] for x, y, z in zip(df_DP["Peptide"], df_DP["pos_peptide"], df_DP["aa_dest"])]
    
    # Diagnostic
    # df_DP.to_csv(os.path.join(dir_path_out, "df_DP_test" + dt.now().strftime("_v%Y%m%d") + ".csv"))
    
    #%% Dynamic subs filtering
    # Full section adapted from write_skyline_input_file.py    
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
        protein_filter = 'Gene.str.match("' + str(args.protein) + '")'
        list_of_filters.append(protein_filter)
        print(INFO_TEXT +
            'FILTERING ... Filter in peptides from ' + str(args.protein) + '.' +
            ENDC_TEXT)
    
    if args.subs_in:
        subs_in = args.subs_in
        subs_in_listoflists = [list(x) for x in subs_in.split("+")]
        subs_in_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_in_listoflists]
        subs_in_filter = 'Sub in @subs_in_list'
        list_of_filters.append(subs_in_filter)
        print(INFO_TEXT +
              'FILTERING ... Filter in substitutions: ' + ', '.join(subs_in_list) +
              ENDC_TEXT)
        
    if args.subs_out:
        subs_out = args.subs_out
        subs_out_listoflists = [list(x) for x in subs_out.split("+")]
        subs_out_list = [" to ".join(x) if len(x) == 2 else (x[0]+" to I/L") for x in subs_out_listoflists]
        subs_out_filter = 'Sub not in @subs_out_list'
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
        df_DP.query(query_cond, inplace=True, engine="python")
    
    if len(df_DP) == 0:
        print(ERROR_TEXT + 
              "ERROR ... No entries in subs after filtering, program was exited. Conflicting filters may have been applied." +
              ENDC_TEXT)
        sys.exit()    
    else:
        print(INFO_TEXT +
              "INFO ...", len(df_DP), "entries in subs after filtering." +
              ENDC_TEXT)
    
    # Diagnostic
    df_DP.to_csv(os.path.join(dir_path_out, "df_DP_filtered" + dt.now().strftime("_v%Y%m%d") + ".csv"))
    
    
    #%% Write fasta file
    
    print(INFO_TEXT +
              "INFO ... Prepare fasta file." +
              ENDC_TEXT)
    
    # Remove replicate entries
    
    # Columns: 
    # "Gene"                = protein name
    # "Peptide"             = BP sequence
    # "modified_sequence"   = DP sequence
    # "pos"                 = position of aa sub in the protein
    # "Sub"                 = "X to Y" style info on subs type (I/L format)
    columns_gb = ["Gene", "Peptide", "modified_sequence", "pos", "Sub"] # TODO no codon information!!!!
    subs_gb = df_DP[columns_gb + ["File_raw"]].groupby(columns_gb).count().reset_index() # raw file column only to allow counting
    
    # Rename columns to allow unpivoting
    subs_gb.rename({"modified_sequence": "Seq_DP", "Peptide": "Seq_BP"}, axis=1, inplace=True)
    subs_gb["id"] = subs_gb.index
    # Unpivot table to generate separate lines for BP & DP
    subs_long = pd.wide_to_long(subs_gb, stubnames="Seq_", i="id", j="Peptide_type", suffix=".+").reset_index(0, drop=True).reset_index(drop=False)
    # Sort unpivoted table in intuitive order (protein > sub position > sub type > BP vs DP)
    subs_long.sort_values(by = ["Gene", "pos", "Sub", "Peptide_type"], inplace=True)
    # Drop duplicate entries
    subs_long.drop_duplicates(subset=["Gene", "Peptide_type", "Seq_"], inplace=True, ignore_index=True)
    
    # Create number for each unique BP and fill down to mark BP & corresponding DP(s)
    subs_long["Num"] = np.nan # Creates peptide number for easier orientation ins Skyline
    subs_long.loc[subs_long["Peptide_type"] == "BP", "Num"] = subs_long[subs_long["Peptide_type"] == "BP"].reset_index().index + 1
    subs_long["Num"] = subs_long["Num"].ffill()
    
    # Write fasta entries per row
    subs_long["fasta"] = subs_long.apply(lambda row: fasta_line(row), axis=1)
    
    subs_long.to_csv(os.path.join(dir_path_out, "subs_long" + dt.now().strftime("_v%Y%m%d") + ".csv"))
    
    # Write fasta file
    print(INFO_TEXT +
              "INFO ... Write fasta file." +
              ENDC_TEXT)
    
    filename1 = dir_path_out + "fasta_skyline" + dt.now().strftime("_v%Y-%m-%d_%H-%M-%S") + ".fa"
    
    f = open(filename1, "a+")
    for index, row in subs_long.iterrows():
        f.write(row["fasta"])
    f.close()
    
    print(OKGREEN_TEXT + 
          "INFO ... Exported fasta file successfully.",
          ENDC_TEXT)
    
    #%% Write ssl file from original df (psm sum file)
    
    print(INFO_TEXT + 
          "INFO ... Prepare ssl file.",
          ENDC_TEXT)
    
    output_DP = pd.DataFrame({
        "file": df_DP["File_raw"],
        "scan": df_DP["Scan"],
        "charge": df_DP["Charge"],
        "sequence": df_DP["modified_sequence"],
        "Score_type" : "Percolator qvalue", #TODO
        "score": df_DP["Qvalue"] # TODO which score
        })
    
    
    # Create list of unique BP sequence that match to identified DPs
    BP_list = list(df_DP["Peptide"].unique())
    # filter DP out of psm to avoid misassigned hits
    df_BP = df.loc[[len(x) == 0 for x in df["mods_clean2"]]].copy() # filter for BP only (no items in mods_clean2 list)    
    # Filter psm for BP
    df_BPonly = df_BP[df_BP["Peptide"].isin(BP_list)]
    
    output_BP = pd.DataFrame({
        "file": df_BPonly["File_raw"],
        "scan": df_BPonly["Scan"],
        "charge": df_BPonly["Charge"],
        "sequence": df_BPonly["Peptide"],
        "Score_type" : "Percolator qvalue", #TODO
        "score": df_BPonly["Qvalue"] #TODO
        })
    
    # this was in write_skyline_input_file.py
    # maybe empty charge values created problems?
    # could be re-evaluated and potentially removed
    output_BP = output_BP[output_BP["charge"] != 0] 
    
    # Prepare final combined table for ssl export    
    output = pd.concat([output_BP, output_DP], ignore_index=True)
    output["sequence"] = output["sequence"].str.replace("C", "C[+57.0]")
    
    # Write ssl file
    print(INFO_TEXT +
              "INFO ... Write ssl file." +
              ENDC_TEXT)
    
    filename2 = dir_path_out + "skyline_input" + dt.now().strftime("_v%Y-%m-%d_%H-%M-%S") + ".ssl"
    output.to_csv(filename2, sep="\t", index=False)
    
    print(OKGREEN_TEXT + 
          "INFO ... Exported ssl file successfully.",
          ENDC_TEXT)
    
    print(OKGREEN_TEXT + 
          "INFO ... Script finished. Good luck",
          ENDC_TEXT)
    
#%%

if __name__ == '__main__':
    main()