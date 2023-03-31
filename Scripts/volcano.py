"""
Generate volcano plots from wild type proteomic changes
"""
import argparse
import warnings
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from sanbomics.plots import volcano

mpl.use('Agg')
warnings.simplefilter('ignore', UserWarning)  # Ignore the warning about plotting with a non-GUI backend


def main():
    parser = argparse.ArgumentParser(description="Generate and save volcano plots. The output filenames are derived "
                                                 "from the input data filename. Generates one plot for 5 min changes "
                                                 "and one for 60 min changes.")
    parser.add_argument("--data", type=str, required=True,
                        help="The file with the q-values and log2 fold changes with the proteomic data.")
    # These should match the thresholds used when generating the data tables
    parser.add_argument("--qvalThresh", type=float, default=0.1,
                        help="The q-value cutoff to use for significance.")
    parser.add_argument("--foldThresh", type=float, default=1.5,
                        help="The fold change cutoff to use for significance.")

    # Get command line args
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    rows = len(df)
    print(f'Read {rows} rows')

    # The Boolean significance calls are not used so drop the column
    df = df.drop(columns='Significant')

    # Iterate over time points
    for time in ['5', '60']:
        # Select the first four columns and the 5 min columns
        pattern = time + ' min'
        time_df = df.iloc[:, 0:4].join(df.loc[:, df.columns.str.contains(pat=pattern)])

        # Select exp 1 and 2 data and rename to remove experiment labels
        exp1_df = time_df.iloc[:, 0:4].join(time_df.loc[:, time_df.columns.str.contains(pat='Exp1')])
        exp1_df.columns = exp1_df.columns.str.replace(f'Exp1 {time} min ', '')

        exp2_df = time_df.iloc[:, 0:4].join(time_df.loc[:, time_df.columns.str.contains(pat='Exp2')])
        exp2_df.columns = exp2_df.columns.str.replace(f'Exp2 {time} min ', '')

        # Reshape the data frame so that the values from both replicates are stacked vertically
        combined_df = pd.concat([exp1_df, exp2_df], ignore_index=True)
        assert combined_df.shape == (rows * 2, 6)

        # Create the output filename by adding the time point to the input data filename
        # The .png and .svg file extension is added automatically
        inpath = Path(args.data)
        outpath = str(Path(inpath.parent, inpath.stem)) + time + 'min'

        # Use to_label=[] to disable labels
        volcano(combined_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[],
                pval_thresh=args.qvalThresh, log2fc_thresh=np.log2(args.foldThresh), pvalue_label='q-value', alpha=0.5,
                linewidth=0.5, symmetric_x_axis=True, colors=['black', 'lightgrey'], save=outpath, legend=False)


if __name__ == "__main__":
    main()
