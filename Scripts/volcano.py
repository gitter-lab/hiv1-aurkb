"""
Generate volcano plots from wild type proteomic changes
"""

import pandas as pd
from sanbomics.plots import volcano

df = pd.read_csv('../Results/base_analysis/AllProteinData.csv')
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
assert combined_df.shape == (rows*2, 6)

print(combined_df.columns)
print(combined_df.head())
print(combined_df)
combined_df.to_csv('../Results/base_analysis/AllProteinDataCombined.csv')

# TODO read the thresholds as arguments
# Use to_label=[] to disable labels
volcano(combined_df, log2fc='log2 fold change', pvalue='q-value', symbol='Uniprot', to_label=[], pval_thresh=0.1,
        log2fc_thresh=0.58496250072, pvalue_label='q-value', alpha=0.5, linewidth=0.5, symmetric_x_axis=True,
        colors=['black', 'lightgrey'])
