#Arguments are: input_prefix output_prefix foldchange_threshold qvalue_threshold

input_prefix=${1:-../Data/mass_spec/csvs/}         # Input prefix to add to all input files.
output_prefix=${2:-../Results/base_analysis/}   # Output prefix to add to all output files.
foldchange_threshold=${3:-1.5}   # Fold change threshold
qvalue_threshold=${4:-0.1} # q-value threshold

#Format data and create enrichment inputs
python mass_spec_analysis.py ${input_prefix}protein_plex1_raw.csv ${input_prefix}protein_plex3_raw.csv ${input_prefix}protein_plex2_raw.csv  ${input_prefix}phospho_plex1_raw.csv ${input_prefix}phospho_plex3_raw.csv ${input_prefix}phospho_plex2_raw.csv --verbose --outName=${output_prefix} --foldThresh=${foldchange_threshold} --qvalThresh=${qvalue_threshold};

#Run limma
for i in `ls -d1 $PWD/${output_prefix}exp*`; do echo $i; Rscript limma_test.R $i; done;

#Create prize lists from limma results
python mass_spec_analysis.py ${input_prefix}protein_plex1_raw.csv ${input_prefix}protein_plex3_raw.csv ${input_prefix}protein_plex2_raw.csv  ${input_prefix}phospho_plex1_raw.csv ${input_prefix}phospho_plex3_raw.csv ${input_prefix}phospho_plex2_raw.csv --verbose --outName=${output_prefix} --foldThresh=${foldchange_threshold} --qvalThresh=${qvalue_threshold};


