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
                                                 "from the input data filename. Generates plots for the 5 min changes "
                                                 "and the 60 min changes. At each time point the individual "
                                                 "experiments are plotted separately and also combined.")
    parser.add_argument("--data", type=str, required=True,
                        help="The file with the q-values and log2 fold changes for the proteomic data.")
    # These should match the thresholds used when generating the data tables
    parser.add_argument("--qvalThresh", type=float, default=0.1,
                        help="The q-value cutoff to use for significance.")
    parser.add_argument("--foldThresh", type=float, default=1.5,
                        help="The fold change cutoff to use for significance. Not yet on log2 scale.")

    # Get command line args
    args = parser.parse_args()

    log2_thresh = np.log2(args.foldThresh)

    # Create the output filename by adding the time point to the input data filename
    # The .png and .svg file extension is added automatically
    inpath = Path(args.data)
    outpath_base = str(Path(inpath.parent, inpath.stem))

    df = pd.read_csv(args.data)
    rows = len(df)
    print(f'Read {rows} rows')

    # The Boolean significance calls are not used so drop the column
    df = df.drop(columns='Significant')

    # Get the maximum -log10 q-value over all time points and experiments to set the y-axis limits
    q_df = df.loc[:, df.columns.str.contains(pat='q-value')]
    assert q_df.shape[1] == 4
    max_q = -np.log10(q_df.min().min())
    print(f'Maximum -log10 q-value: {max_q}')
    max_q = max_q * 1.025  # add extra padding

    # Get the maximum magnitude -log10 fold change over all time points and experiments to set the x-axis limits
    fc_df = df.loc[:, df.columns.str.contains(pat='log2 fold change')]
    assert fc_df.shape[1] == 4
    max_fc = fc_df.abs().max().max()
    print(f'Maximum magnitude log2 fold change: {max_fc}')
    max_fc = max_fc * 1.05  # add extra padding

    # Iterate over time points
    for time in ['5', '60']:
        # Select the first four columns and the 5 min columns
        pattern = time + ' min'
        time_df = df.iloc[:, 0:4].join(df.loc[:, df.columns.str.contains(pat=pattern)])

        # Select exp 1 and 2 data and rename to remove experiment labels
        exp1_df = time_df.iloc[:, 0:4].join(time_df.loc[:, time_df.columns.str.contains(pat='Exp1')])
        exp1_df.columns = exp1_df.columns.str.replace(f'Exp1 {time} min ', '')

        outpath1 = outpath_base + time + 'minOnlyExp1'
        volcano(exp1_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[],
                pval_thresh=args.qvalThresh, log2fc_thresh=log2_thresh, pvalue_label='q-value', alpha=0.5,
                linewidth=0.5, symmetric_x_axis=True, colors=['black', 'lightgrey'], save=outpath1, legend=False)

        exp2_df = time_df.iloc[:, 0:4].join(time_df.loc[:, time_df.columns.str.contains(pat='Exp2')])
        exp2_df.columns = exp2_df.columns.str.replace(f'Exp2 {time} min ', '')

        outpath2 = outpath_base + time + 'minOnlyExp2'
        volcano(exp2_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[],
                pval_thresh=args.qvalThresh, log2fc_thresh=log2_thresh, pvalue_label='q-value', alpha=0.5,
                linewidth=0.5, symmetric_x_axis=True, colors=['black', 'lightgrey'], save=outpath2, legend=False)

        # Reshape the data frame so that the values from both replicates are stacked vertically
        combined_df = pd.concat([exp1_df, exp2_df], ignore_index=True)
        assert combined_df.shape == (rows * 2, 6)

        outpath_combined = outpath_base + time + 'min'
        # Use to_label=[] to disable labels
        # Only use the manual x and y limits when combining data from both experiments
        volcano(combined_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[],
                pval_thresh=args.qvalThresh, log2fc_thresh=log2_thresh, pvalue_label='q-value', alpha=0.5,
                linewidth=0.5, symmetric_x_axis=True, x_max=max_fc, y_max=max_q, colors=['black', 'lightgrey'],
                save=outpath_combined, legend=False)


if __name__ == "__main__":
    main()
