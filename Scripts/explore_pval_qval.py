import seaborn as sns
import pandas as pd
import numpy as np

results_dir = '../Results/base_analysis/'

df = pd.read_csv(results_dir + 'expPhos2_05.csv_tmpLimma')

# -log10 transform the p-values and q-values
df['neg_log10_pval'] = -np.log10(df['pval'])
df['neg_log10_qval'] = -np.log10(df['qval'])

# histograms of original and -log10 p-values and q-values
fig = sns.histplot(data=df, x='pval').get_figure()
fig.savefig(results_dir + 'expPhos2_05_pval_hist.png')
fig.clf()
fig = sns.histplot(data=df, x='qval').get_figure()
fig.savefig(results_dir + 'expPhos2_05_qval_hist.png')
fig.clf()
fig = sns.histplot(data=df, x='neg_log10_pval').get_figure()
fig.savefig(results_dir + 'expPhos2_05_neg_log10_pval_hist.png')
fig.clf()
fig = sns.histplot(data=df, x='neg_log10_qval').get_figure()
fig.savefig(results_dir + 'expPhos2_05_neg_log10_qval_hist.png')
fig.clf()

# scatterplot of -log10 p-values vs q-values
fig = sns.scatterplot(df, x='neg_log10_pval', y='neg_log10_qval').get_figure()
fig.savefig(results_dir + 'expPhos2_05_neg_log10_qval_vs_pval.png')
fig.clf()
