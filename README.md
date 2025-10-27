# skyline_pipeline

This script enables the quantitative analysis of amino acid subtitutions from mass spectrometry data in Skyline by writing files for importing peptide search results (.ssl and .fasta files). It relies on an initial analysis of mass spectrometry raw files using the MaxQuant Dependent Peptide search function and the subsequent filtering performed by the collection of scripts published by [Mordret et al. (2019)](https://doi.org/10.1016/j.molcel.2019.06.041).

## Available scripts
`write_skyline_input_file.py` is the default script and is operated from the terminal. 
`write_skyline_input_file_addon.py` takes entries from a curated list of peptides that can generally not be processed with the Dependent Peptide feature but could be identified in the MaxQuant search (e.g. peptides derived from missense errors leading to altered trypsin cleavage patterns). This script does not support loggers or operation from the terminal.
`write_skyline_input_file_manual_input.py` allows filtering of the missense peptide list (subs file, derived from Mordret script) prior to generating the Skyline input files. Deprecated, this function is now integrated into the main script `write_skyline_input_file.py`.

## Example 
Create Skyline input files with only near-cognate substitutions in the protein EF-Tu that are not found in the unimod database (with the exception off D to E).
```
python write_skyline_input_file.py -m -p "tuf[AB]" -ft "(danger == False or substitution == 'D to E')" <path to input dir>
```

For information on all filtering options, run
```
python write_skyline_input_file.py -h
```