# skyline_pipeline

This script enables the quantitative analysis of amino acid subtitutions from mass spectrometry data in Skyline by writing files for importing peptide search results (.ssl and .fasta files). It relies on an initial analysis of mass spectrometry raw files using one of the following two search engines: 
1. MaxQuant: Modified peptides are identified using dependent peptide searching. Subsequent filtering for amino acid substitutions is required and needs to be performed using the collection of scripts published by [Mordret et al. (2019)](https://doi.org/10.1016/j.molcel.2019.06.041), using either the [original script collection](https://github.com/ernestmordret/substitutions) or an [adapted version](https://github.com/nfreyer/substitutions) tailored for the present workflow.
2. MSFragger: Modified peptides are identified using an unrestricted mass offset search. The FragPipe workflow file and corresponding mass offset list are provided here.

## Available scripts
`write_maxquant_to_skyline.py` uses the results of the MaxQuant dependent peptide search and substitution filtering. Operated from the terminal.

`write_fragpipe_to_skyline.py` uses the results of the mass offset search in MSFragger and performs filtering for substitutions itself. Operated from the terminal.

`write_skyline_input_file_addon.py` takes entries from a curated list of peptides that can generally not be processed with the Dependent Peptide feature but could be identified in the MaxQuant search (e.g. peptides derived from missense errors leading to altered trypsin cleavage patterns). This script does not support loggers or operation from the terminal.

## Example 
Create Skyline input files with only near-cognate substitutions in the protein EF-Tu that are not found in the unimod database (with the exception off D to E).
```
python write_maxquant_to_skyline.py -m -p "tuf[AB]" -ft "(danger == False or substitution == 'D to E')" <path to input dir>
```

For information on all filtering options, run
```
python write_maxquant_to_skyline.py -h
```
or 
```
python write_fragpipe_to_skyline.py -h
```