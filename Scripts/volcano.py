"""
Generate volcano plots from wild type proteomic changes
"""
import argparse

import numpy as np
import pandas as pd
from sanbomics.plots import volcano


def main():
    parser = argparse.ArgumentParser(description="Generate and save volcano plots. The output filenames are derived "
                                                 "from the input data filename. Generates one plot for 5 min changes"
                                                 "and one for 60 min changes.")
    parser.add_argument("--data", type=str, required=True,
                        help="This is the file with the q-values and log2 fold changes with the proteomic data.")
    # These should match the thresholds used when generating the data tables
    parser.add_argument("--qvalThresh", type=float, default=0.1,
                        help="This is the q-value cutoff to use for significance.")
    parser.add_argument("--foldThresh", type=float, default=1.5,
                        help="This is the fold change cutoff to use for significance.")

    # Get command line args
    args = parser.parse_args()

    df = pd.read_csv(args.data)
    rows = len(df)
    print('Read {rows} rows')

    df = df.drop(columns='Significant')

    # TODO iterate over time points
    time = '5 min'

    # Select the first four columns and the 5 min columns
    df5 = df.iloc[:, 0:4].join(df.loc[:, df.columns.str.contains(pat=time)])

    # Select exp 1 and 2 data and rename to remove experiment labels
    exp1_df = df5.iloc[:, 0:4].join(df5.loc[:, df5.columns.str.contains(pat='Exp1')])
    exp1_df.columns = exp1_df.columns.str.replace('Exp1 5 min ', '')

    exp2_df = df5.iloc[:, 0:4].join(df5.loc[:, df5.columns.str.contains(pat='Exp2')])
    exp2_df.columns = exp2_df.columns.str.replace('Exp2 5 min ', '')

    combined_df = pd.concat([exp1_df, exp2_df], ignore_index=True)
    assert combined_df.shape == (rows * 2, 6)

    print(combined_df.columns)
    print(combined_df.head())
    print(combined_df)
    combined_df.to_csv('../Results/base_analysis/AllProteinDataCombined.csv')

    # Use to_label=[] to disable labels
    volcano(combined_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[], pval_thresh=args.qvalThresh,
            log2fc_thresh=np.log2(args.foldThresh), pvalue_label='q-value', alpha=0.5, linewidth=0.5, symmetric_x_axis=True,
            colors=['black', 'lightgrey'])


if __name__ == "__main__":
    main()
