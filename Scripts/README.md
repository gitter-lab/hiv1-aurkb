# Scripts
There are two primary scripts to regenerate data tables and figures
- `run_limma.sh` - the main script that compares time points in the wild type conditions, which was used for nearly all analyses in the manuscript
- `run_limma_control.sh` - a script analyzing -Env and -Genome control conditions

Both scripts follow a similar pattern: they run a Python script to collect values from the `.csv` files, they run limma for differential abundance or phosphorylation testing, and they run the Python script again for post-processing.
The scripts must be run in two stages, however, with the unused half of the script commented out.

All code should be run in the conda environment created by `environment.yml`.

## Other files
- `compare_control_conditions.py` - Python script used only in the control analysis
- `environment.yml` - conda environment file specifying both Python and R dependencies
- `limma_test.R` - R script used to run limma in both analyses
- `mass_spec_analysis.py` - Python script used only in the wild type analysis
- `volcano.py` - Python script for generating figures not currently included in the main `.sh` scripts
