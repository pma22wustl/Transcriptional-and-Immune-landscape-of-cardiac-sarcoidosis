import palantir
from sklearn.preprocessing import StandardScaler

import scanpy as sc
import numpy as np
import pandas as pd

# Plotting 
import matplotlib as mt
import matplotlib.pyplot as plt
import seaborn as sns
import gam

import os
os.environ['R_HOME'] = '/Library/Frameworks/R.framework/Versions/4.0/Resources' #path to your R installation
os.environ['R_USER'] = '/Users/lululai/opt/anaconda3/Lib/site-packages/rpy2' #path depends on where you installed Python.

import rpy2.robjects.packages as rpackages
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)
# Import packages
from rpy2.robjects.packages import importr
base, bnlearn = importr('base'), importr('bnlearn')

# Reset random seed
np.random.seed(5)

palantir_dir = '/Users/lululai/Box Sync/Sarcoidosis Project/Sarc_only_integration/cleaned/Myeloid.new.ident/Palantir/'

counts = pd.read_csv('/Users/lululai/Box Sync/Sarcoidosis Project/Sarc_only_integration/cleaned/Myeloid.new.ident/Palantir/myeloid_RNA_data.txt', sep=',', index_col=0).transpose()
counts = counts.sort_index()

meta_data = pd.read_csv('/Users/lululai/Box Sync/Sarcoidosis Project/Sarc_only_integration/cleaned/Myeloid.new.ident/Palantir/myeloid_meta.csv', index_col=0)
meta_data = meta_data.sort_index()

pca_projections, _ = palantir.utils.run_pca(counts, use_hvg=False)
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res)
tsne = palantir.utils.run_tsne(ms_data)

fig, ax = palantir.plot.plot_tsne(tsne)
plt.savefig(palantir_dir+'tsne.png', dpi=1200)

fig, ax = palantir.plot.plot_tsne_by_cell_sizes(counts, tsne)
plt.savefig('tsne_counts.png', dpi=1200)

imp_df = palantir.utils.run_magic_imputation(counts, dm_res)
palantir.plot.plot_diffusion_components(tsne, dm_res)
plt.savefig('tsne_diff_components.png', dpi=1200)

tsne_with_index = tsne.set_index([meta_data.index])
palantir.plot.plot_cell_clusters(tsne_with_index, meta_data["orig.ident"])
plt.savefig('tsne_samples.png', dpi=1200)

palantir.plot.plot_cell_clusters(tsne_with_index, meta_data["new.ident"])
plt.savefig('tsne_new.ident.png', dpi=1200)

meta_data_mono = meta_data[meta_data["new.ident"] == 6]
meta_data_mono

start_cell = "TWCM-574-574_ATCTACTCACAAGACG"
palantir.plot.highlight_cells_on_tsne(tsne_with_index, [start_cell])
plt.savefig('start_cell.png', dpi=1200)

#terminal_states = pd.Series(['0', '1', '2'], 
#                           index=['Run5_131097901611291', 'Run5_134936662236454', 'Run4_200562869397916'])

# counts = palantir.io.from_csv(palantir_dir+'raw_counts.csv')
# counts = counts.transpose()
# fig, ax = palantir.plot.plot_molecules_per_cell_and_gene(counts)
# fig = plt.gcf()
# fig.savefig(palantir_dir+'pseudotime.png', dpi=1200)
#terminal_states = pd.Series(['resident'], 
#                           index=['M_JC-L3-L3_TAGGCATCATTCCTGC'])
msdata_with_index = ms_data.set_index([meta_data.index])
pr_res = palantir.core.run_palantir(msdata_with_index, start_cell, num_waypoints=500)

pr_res.branch_probs.columns
# lst = []
# for i in list(pr_res.branch_probs.columns):
#     print(meta_data["new.ident"][i])
#     lst.append(meta_data["new.ident"][i])

# pr_res.branch_probs.columns = lst
# pr_res.branch_probs = pr_res.branch_probs.loc[:, lst]

#fdl = harmony.plot.force_directed_layout(dm_res['kernel'], ad.obs.new_ident)

palantir.plot.plot_palantir_results(pr_res, tsne_with_index)
plt.savefig('palantir_summary.png', dpi=1200)

# palantir.plot.plot_gene_expression(imp_df, tsne_with_index, ['TPRG1','GPNMB','LILRB1'])
# plt.savefig('gene_exp_result.png', dpi=1200)

# genes=['TPRG1','GPNMB','LILRB1']
# gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, genes])

#create meta data
df = pd.DataFrame()
df = pr_res.branch_probs.copy()
df['pseudotime'] = np.array(pr_res.pseudotime)
df['entropy'] = np.array(pr_res.entropy)
df['ClusterName'] = list(meta_data["new.ident"])
df['Samples'] = list(meta_data["orig.ident"])
#df['Group'] = list(meta_data.group)
df['ClusterName'] = df['ClusterName'].astype("category")
df['Samples'] = df['Samples'].astype("category")

#saving data
pca_projections.to_csv('myeloid_pca_projections.csv', index=True)
ms_data.to_csv('myeloid_ms_data.csv', index=True)
tsne.to_csv('myeloid_tsne.csv', index=True)
df.to_csv('myeloid_palantir_meta_data.csv', index=True)

#pseudotime vs entropy by cluster names
ax = sns.scatterplot(x="pseudotime", y="entropy", hue="ClusterName", data=df, s=10,linewidth=0, legend="auto", palette="Set2")
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('pseudotime_entropy_annotation.png', bbox_inches='tight', dpi=1200)
plt.clf()

#pseudotime vs entropy by smaples
ax = sns.scatterplot(x="pseudotime", y="entropy", hue="Samples", data=df, s=10,linewidth=0, legend="auto", palette="Set2")
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
plt.savefig('pseudotime_entropy_samples.png', bbox_inches='tight', dpi=1200)
plt.clf()

#pseudotime vs entropy individual sample plots
df_574 = df[df["Samples"]=="TWCM-574-574"]
df_676 = df[df["Samples"]=="TWCM-676-676"]
df_707 = df[df["Samples"]=="TWCM-707-707"]
df_761 = df[df["Samples"]=="TWCM-761-761"]

a_dict={}
#sample_df = np.array(['df_574',df_574],['df_676',df_676], ['df_707',df_707], ['df_761',df_761])
for variable in ["df_574", "df_676", "df_707","df_761"]:
    a_dict[variable] = eval(variable)
mylist = list(a_dict.values())

for i, value in enumerate(mylist):   
    ax = sns.scatterplot(x="pseudotime", y="entropy", hue="Samples", data=mylist[i], s=10,linewidth=0, legend="auto", palette="Set2")
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
    plt.savefig('pseudotime_entropy_samples%s.png'%i, bbox_inches='tight', dpi=1200)
    plt.clf()
    
# sample_df = np.array[['df_574',df_574],['df_676',df_676], ['df_707',df_707], ['df_761',df_761]]

# for i, value in enumerate(sample_df):    
#     ax = sns.scatterplot(x="pseudotime", y="entropy", hue="Samples", data=sample_df[i], s=10,linewidth=0, legend="auto", palette="Set2")
#     #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
#     plt.savefig('pseudotime_entropy_samples%s.png'%i, bbox_inches='tight', dpi=1200)
#     plt.clf()

#gene expression
gene_list = ['HLA-DRA', 'TPRG1', 'GPNMB', 'CHIT1', 'CHI3L1',\
                 'CD44', 'RUNX2', 'PLAUR', 'MMP9', 'STAT1', 'FBP1',\
                     'LYZ', 'VIM', 'CXCL9', 'S100A6','S100A11','ITGAX',\
                         'CCL4L2']
gene_list2 = ["LAMP3", "LAMP2", "CD63", \
             "CD58","CD164", "CD300LF", "CD46", "CD48", "CD52"]
palantir.plot.plot_gene_expression(imp_df, tsne_with_index, gene_list)
plt.savefig('gene_exp.png', bbox_inches='tight', dpi=1200)
plt.clf()

#lineage trends determination
gene_trends3 = palantir.presults.compute_gene_trends(pr_res, imp_df[gene_list], 
                                                     lineages=['TWCM-676-676_AGGGATGGTTACTGAC', 'TWCM-707-707_CCATGTCCAAGCGTAG',
       'TWCM-761-761_ACAGCTAAGGTGATTA', 'TWCM-761-761_TGACAACCATATACGC'])
gene_trends4 = palantir.presults.compute_gene_trends(pr_res, imp_df[gene_list2], lineages=['TWCM-676-676_AGGGATGGTTACTGAC', 'TWCM-707-707_CCATGTCCAAGCGTAG',
       'TWCM-761-761_ACAGCTAAGGTGATTA', 'TWCM-761-761_TGACAACCATATACGC'])

#plot trends
def plot_gene_trends(gene_trends, genes=None):
    """ Plot the gene trends: each gene is plotted in a different panel
    :param: gene_trends: Results of the compute_marker_trends function
    """

    # Branches and genes
    branches = list(gene_trends.keys())
    colors = pd.Series(
        sns.color_palette("Set2", len(branches)).as_hex(), index=branches
    )
    if genes is None:
        genes = gene_trends[branches[0]]["trends"].index

    # Set up figure
    fig = plt.figure(figsize=[7, 3 * len(genes)])
    for i, gene in enumerate(genes):
        ax = fig.add_subplot(len(genes), 1, i + 1)
        for branch in branches:
            trends = gene_trends[branch]["trends"]
            stds = gene_trends[branch]["std"]
            ax.plot(
                trends.columns, trends.loc[gene, :], color=colors[branch], label=branch
            )
            ax.set_xticks([0, 1])
            ax.fill_between(
                trends.columns,
                trends.loc[gene, :] - stds.loc[gene, :],
                trends.loc[gene, :] + stds.loc[gene, :],
                alpha=0.1,
                color=colors[branch],
            )
            ax.set_title(gene)
        # Add legend
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    sns.despine()


plot_gene_trends(gene_trends3)
plt.savefig('lineage_trends_1.png', bbox_inches='tight', dpi=1200)
plt.clf()

plot_gene_trends(gene_trends4)
plt.savefig('lineage_trends_2.png', bbox_inches='tight', dpi=1200)
plt.clf()

#heatmaps
def plot_gene_trend_heatmaps(gene_trends):
    """ Plot the gene trends on heatmap: a heatmap is generated or each branch
    :param: gene_trends: Results of the compute_marker_trends function
    """

    # Plot height
    branches = list(gene_trends.keys())
    genes = gene_trends[branches[0]]["trends"].index
    height = 0.7 * len(genes) * len(branches)

    #  Set up plot
    fig = plt.figure(figsize=[14, height/2])
    #fig, ax = plt.subplots(2,2, figsize=[28, 3 * len(genes)])
    for i, branch in enumerate(branches):
        ax = fig.add_subplot(2, 2, i + 1)

        # Standardize the matrix
        mat = gene_trends[branch]["trends"]
        mat = pd.DataFrame(
            StandardScaler().fit_transform(mat.T).T,
            index=mat.index,
            columns=mat.columns,
        )
        sns.heatmap(mat, xticklabels=False, ax=ax, cmap=mt.cm.Spectral_r)
        ax.set_title(branch, fontsize=12)
        
plot_gene_trend_heatmaps(gene_trends3)
plt.savefig('heatmap1.png', dpi=300)
plt.clf()

plot_gene_trend_heatmaps(gene_trends4)
plt.savefig('heatmap2.png', dpi=300)
plt.clf()
