# skyline_pipeline

This script enables the quantitative analysis of amino acid subtitutions from mass spectrometry data in Skyline by writing files for importing peptide search results (.ssl and .fasta files). It relies on an initial analysis of mass spectrometry raw files using the MaxQuant Dependent Peptide search function and the subsequent filtering performed by the collection of scripts published by Mordret et al. (2019).

## Example: 
Create Skyline input files with only near-cognate substitutions in the protein EF-Tu that are not found in the unimod database (with the exception off D to E).
	python write_skyline_input_file.py -m -p "tuf[AB]" -ft "(danger == False or substitution == 'D to E')" <path to input dir>

For information on all filtering options, run
	python write_skyline_input_file.py -h