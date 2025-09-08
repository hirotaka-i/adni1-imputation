# Input: plink1.9 pca results of the study merged with the 1000 genome project
# Output1: Ancestry labeling determined by +/- 6SD of PC1-5 of the reference panel
# Output2: PCs with the 1000 genome project
# Output3: Plots for the whole cohort and the European cohort. 

# PCA results visualization and population inference
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import argparse
import os
import numpy as np
from scipy.spatial import distance
import matplotlib.cm as cm

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description="calculate pca from plink file")

# Define the command-line arguments
parser.add_argument("--pc-prefix", type=str, help="prefix of the plink pca output files", required=True)
parser.add_argument("--ref-label", type=str, help="population label file for the reference panel", default="all_hg38_filtered_chrpos_pop.txt")
parser.add_argument("--ref-label-col", type=str, help="population label column for the reference panel", required=True)
parser.add_argument("--out-prefix", type=str, help="prefix of the output files", default="genetic_ancestry")
parser.add_argument("--split-method", type=str, choices=['none', 'mahalanobis', 'sd'], default='none',
                    help="Method for ancestry split: none, mahalanobis, or sd (mean+/-4SD)")
parser.add_argument("--eur-aj-sep", action='store_true', default=False,
                    help="If set, split EUR into EUR-AJ and EUR-nonAJ using eur_pop_dict")


# Parse the command-line arguments
args = parser.parse_args()
pc_prefix = args.pc_prefix
ref_label = args.ref_label
ref_label_col = args.ref_label_col
out_prefix = args.out_prefix
split_method = args.split_method
eur_aj_sep = args.eur_aj_sep

# # Example values for testing
# pc_prefix = 'temp/EUR/merge_ref_proj_out/study_vs_ref.combined'
# ref_label = 'data/all_hg38_filtered_chrpos_pop.txt'
# out_prefix = ''
# ref_label_col = 'Population'
# study_list = '_EUR.list'
# ref_filter = 'EUR'
# split_ancestry = True

# make sure output directory exists
output_dir = os.path.dirname(out_prefix)
if output_dir and not os.path.exists(output_dir):
    os.makedirs(output_dir) 

# Population dictionary for grouping
pop_dict = {
    'AMR': ['MXL','CLM','PEL','PUR'],
    'EAS': ['JPT','CDX','CHB','CHS','KHV','CHD'],
    'EUR': ['TSI','IBS','GBR','CEU', 'AJ', 'FIN'],
    'SAS': ['PJL','ITU','STU','GIH','BEB'],
    'AFR': ['GWD','MSL','ESN','GWJ','YRI','LWK','GWF','GWW'],
    'AAC': ['ASW','ACB'],
}
eur_pop_dict = {
    'EUR-nonAJ': ['TSI','IBS','GBR','CEU', 'FIN'],
    'EUR-AJ': ['AJ'],
}

# Flexible color assignment for populations
def get_colors(n, cmap_name='tab20'):
    cmap = plt.colormaps[cmap_name]
    return [cmap(i) for i in np.linspace(0, 1, n)]

# Population figure creation
def make_and_save_population_figure(df, label_col, out_prefix, title, x_pc='PC1', y_pc='PC2', cmap_name='tab20'):
    labs = df[label_col].unique()
    colors = get_colors(len(labs), cmap_name)
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    legend_lines = []
    for i, (j, group) in enumerate(df.groupby(label_col)):
        if j == 'Study':
            ax.scatter(x=x_pc, y=y_pc, color=colors[i], label=j, s=10, alpha=1, data=group)
            n_study = group.shape[0]
        else:
            sns.kdeplot(x=x_pc, y=y_pc, n_levels=4, ax=ax, color=colors[i], alpha=0.8, data=group)
            legend_lines.append(mlines.Line2D([], [], color=colors[i], label=j))
    ax.legend(handles=legend_lines, title='Label', loc='upper right')
    plt.title(f'{title},n={n_study}:{x_pc}_vs_{y_pc}')
    fig.savefig(f'{out_prefix}_{title}_{x_pc}x{y_pc}.png', dpi=300, format='png')
    plt.close(fig)

# Ancestry assignment using Mahalanobis distance
def assign_population_mahalanobis(row, pop_stats, threshold=6):
    if row['label_group'] != 'Study':
        return 'REF'
    x = row[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']].values
    min_dist = float('inf')
    min_pop = 'OTHER'
    for pop, stats in pop_stats.items():
        dist = distance.mahalanobis(x, stats['mean'], stats['inv_cov'])
        if dist < min_dist:
            min_dist = dist
            min_pop = pop
    if min_dist < threshold:
        return min_pop
    else:
        return 'OTHER'

# Ancestry assignment using mean+/-4SD
def assign_population_sd(row, pop_stats):
    if row['label_group'] != 'Study':
        return 'REF'
    for pop, stats in pop_stats.items():
        in_range = True
        for i, pc in enumerate(['PC1', 'PC2', 'PC3', 'PC4', 'PC5']):
            val = row[pc]
            mean = stats['mean'][i]
            std = stats['std'][i]
            if not (mean - 4*std < val < mean + 4*std):
                in_range = False
                break
        if in_range:
            return pop
    return 'OTHER'

# Scree plot
dfpcaVal = pd.read_csv(f'{pc_prefix}.eigenval', header=None)
dfpcaVal.plot.line()
plt.title('Scree plot')
plt.savefig(f'{out_prefix}_screeplot.png', dpi=300, format='png')

# PCA plot
# Read the population labeling
t=pd.read_csv(ref_label, sep='\t')
t['label'] = t[ref_label_col].fillna('UNKNOWN')

if eur_aj_sep:
    def eur_group(x):
        for group, ids in eur_pop_dict.items():
            if x in ids:
                return group
        # fallback to original pop_dict
        for pop, ids in pop_dict.items():
            if x in ids:
                return pop
        return 'OTHER'
    t['label_group'] = t['label'].apply(eur_group)
else:
    t['label_group'] = t['label'].apply(lambda x: next((pop for pop, ids in pop_dict.items() if x in ids), 'OTHER'))

# Merge with eigenvec
d=pd.read_csv(f'{pc_prefix}.eigenvec', sep=r'\s+')
df=pd.merge(d, t[['IID', 'label', 'label_group']], on=['IID'], how='left')
df['label'] = df['label'].fillna('Study')
df['label_group']=df['label_group'].fillna('Study')
print(df.groupby(['label_group', 'label']).size())

# General figures
make_and_save_population_figure(
    df=df,
    label_col='label_group',
    out_prefix=out_prefix,
    title='Whole group',
    x_pc='PC1',
    y_pc='PC2'
)

make_and_save_population_figure(
    df=df,
    label_col='label',
    out_prefix=out_prefix,
    title='Whole ancestry',
    x_pc='PC1',
    y_pc='PC2'
)

# Infer ancestry using selected split method
if split_method != 'none':
    pop_stats = {}
    pops_in_ref = df[df.label_group != 'Study']['label_group'].unique()
    for pop in pops_in_ref:
        pop_df = df[df['label_group'] == pop][['PC1', 'PC2', 'PC3', 'PC4', 'PC5']]
        mean_vec = pop_df.mean().values
        cov_mat = np.cov(pop_df.values, rowvar=False)
        std_vec = pop_df.std().values
        # Regularize covariance matrix if singular
        if np.linalg.det(cov_mat) == 0:
            cov_mat += np.eye(cov_mat.shape[0]) * 1e-6
        pop_stats[pop] = {'mean': mean_vec, 'cov': cov_mat, 'inv_cov': np.linalg.inv(cov_mat), 'std': std_vec}

    if split_method == 'mahalanobis':
        df['InfPop'] = df.apply(assign_population_mahalanobis, axis=1, pop_stats=pop_stats)
    elif split_method == 'sd':
        df['InfPop'] = df.apply(assign_population_sd, axis=1, pop_stats=pop_stats)

    print(df.InfPop.value_counts())
    pops_in_study = df[df.label == 'Study']['InfPop'].unique()
    # save df with inferred population
    df.loc[df.InfPop!='REF', ['#FID', 'IID', 'InfPop'] + [f'PC{i+1}' for i in range(10)]].to_csv(f'{out_prefix}_all_pca.csv', index=False)
    # save the pop lists and generate plots
    for continent in pops_in_study:
        t = df.loc[df.InfPop==continent, ['#FID', 'IID']]
        t.to_csv(f'{out_prefix}_{continent}.list', sep='\t', index=False, header=False)
        print(continent, t.shape[0], 'list saved to', f'{out_prefix}_{continent}.list')
        if continent != 'OTHER':
            for y_pc in ['PC2', 'PC3']:
                make_and_save_population_figure(
                    df=df[(df.InfPop==continent) | (df.label_group==continent)],
                    label_col='label',
                    out_prefix=out_prefix,
                    title=f'{continent}',
                    x_pc='PC1',
                    y_pc=y_pc
                )
