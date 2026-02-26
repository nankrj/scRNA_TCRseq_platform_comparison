import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import scipy.io as sio
from scipy.stats import median_abs_deviation
import seaborn as sns
import scanpy.external as sce
import matplotlib.pyplot as plt
import anndata as ad
from matplotlib_venn import venn2


#### Loading the Parse data - called adata
# The DGE_filtered folder contains the expression matrix, genes, and files
adata = sc.read_mtx('/path/to/folder/with/count_matrix.mtx')

adata.write(h5ad_path + 'adata_PARSE.h5ad')
# adata = sc.read(obj_save_path + 'adata_obj1.h5ad')

# reading in gene and cell data
gene_data = pd.read_csv(input_data_path + 'all-sample/DGE_filtered/all_genes.csv')
cell_meta = pd.read_csv(input_data_path + 'all-sample/DGE_filtered/cell_metadata.csv')

# find genes with nan values and filter
gene_data = gene_data[gene_data.gene_name.notnull()]
notNa = gene_data.index
notNa = notNa.to_list()

# remove genes with nan values and assign gene names
adata = adata[:,notNa]
adata.var = gene_data
adata.var.set_index('gene_name', inplace=True)
adata.var.index.name = None
adata.var_names_make_unique()

# add cell meta data to anndata object
adata.obs = cell_meta
adata.obs.set_index('bc_wells', inplace=True)
adata.obs.index.name = None
adata.obs_names_make_unique(join = '_')

# adding patient information
group_info = {
    'ID1_Blood_Tcell_PARSE': 'ID1',
    'ID2_Blood_Tcell_PARSE': 'ID2',
    'ID3_Blood_Tcell_PARSE': 'ID3',
    'ID4_Blood_Tcell_PARSE': 'ID4',
    'ID5_Blood_Tcell_PARSE': 'ID5',
    'ID6_Blood_Tcell_PARSE': 'ID6',
    'ID7_Blood_Tcell_PARSE': 'ID7'
}

figure_label = {
    'ID1_Blood_Tcell_PARSE': '1',
    'ID2_Blood_Tcell_PARSE': '2',
    'ID3_Blood_Tcell_PARSE': '3',
    'ID4_Blood_Tcell_PARSE':  '4',
    'ID5_Blood_Tcell_PARSE': '5',
    'ID6_Blood_Tcell_PARSE':  '6',
    'ID7_Blood_Tcell_PARSE':  '7'
}
adata.obs['patient'] = adata.obs['sample'].map(group_info)
adata.obs['fig_label'] = adata.obs['sample'].map(figure_label)

adata = adata[adata.obs.patient.notna()]


#### Calculating QC metrics

# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)

#Outlier function
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


adata.obs["outlier"] = (
    is_outlier(adata, "log1p_total_counts", 5) |
    is_outlier(adata, "log1p_n_genes_by_counts", 5) |
    is_outlier(adata, "pct_counts_in_top_20_genes", 5) |
    (adata.obs["pct_counts_mt"] > 10)
)
adata.obs.outlier.value_counts()


dt = adata.obs[['fig_label','outlier']].value_counts().reset_index(name = 'count').pivot(index = 'fig_label', columns = 'outlier', values = 'count')
dt2 = dt.div(dt.sum(axis = 1), axis = 0)*100

#### Loading the 10X data - called tadata - calculating QC and finding outliers

input_data_path = 'path/to/folder/with/10X/multi/output/'

filenames = ['ID1_Blood_Tcell_10X', 'ID2_Blood_Tcell_10X', 'ID3_Blood_Tcell_10X', 'ID4_Blood_Tcell_10X', 'ID5_Blood_Tcell_10X',
            'ID6_Blood_Tcell_10X', 'ID7_Blood_Tcell_10X']

tadatas = {}

for filename in filenames:
    sample_adata = sc.read_10x_h5(input_data_path + filename + '/outs/per_sample_outs/' + filename + '/count/sample_feature_bc_matrix.h5')
    sample_adata.var_names_make_unique()
    
    fca = pd.read_csv(input_data_path + filename + '/outs/per_sample_outs/' + filename + '/vdj_t/filtered_contig_annotations.csv')
    fca = fca.drop_duplicates(subset='barcode')
    fca = fca[['barcode', 'raw_clonotype_id']]
    fca.rename(columns={"raw_clonotype_id": "clonotype_id"}, inplace=True)
    
    ct = pd.read_csv(input_data_path + filename + '/outs/per_sample_outs/' + filename + '/vdj_t/clonotypes.csv')
    
    mdf = fca.merge(ct, on = 'clonotype_id', how='inner')
    mdf.set_index('barcode', inplace = True)
    mdf.index.name = None

    obs_df = sample_adata.obs.merge(mdf, left_index = True, right_index = True, how = 'left')

    sample_adata.obs = obs_df
    #sample_adata.obs.set_index('barcode', inplace=True)
    #sample_adata.obs.index.name = None
    sample_adata.obs_names_make_unique()

    # mitochondrial genes
    sample_adata.var["mt"] = sample_adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(sample_adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)


    sample_adata.obs["outlier"] = (
    is_outlier(sample_adata, "log1p_total_counts", 5) |
    is_outlier(sample_adata, "log1p_n_genes_by_counts", 5) |
    is_outlier(sample_adata, "pct_counts_in_top_20_genes", 5) |
    (sample_adata.obs["pct_counts_mt"] > 10)
    )

    tadatas[filename] = sample_adata

tadata = ad.concat(tadatas, label='sample')
tadata.obs_names_make_unique()
print(tadata.obs['sample'].value_counts())

group_info = {
    'ID1_Blood_Tcell_10X': 'ID1',
    'ID2_Blood_Tcell_10X': 'ID2',
    'ID3_Blood_Tcell_10X': 'ID3',
    'ID4_Blood_Tcell_10X': 'ID4',
    'ID5_Blood_Tcell_10X': 'ID5',
    'ID6_Blood_Tcell_10X': 'ID6',
    'ID7_Blood_Tcell_10X': 'ID7'
}


tadata.obs['patient'] = tadata.obs['sample'].map(group_info)

fig_info = {
    'ID1_Blood_Tcell_10X': '1',
    'ID2_Blood_Tcell_10X': '2',
    'ID3_Blood_Tcell_10X': '3',
    'ID4_Blood_Tcell_10X': '4',
    'ID5_Blood_Tcell_10X': '5',
    'ID6_Blood_Tcell_10X': '6',
    'ID7_Blood_Tcell_10X': '7'
}


tadata.obs['fig_label'] = tadata.obs['sample'].map(fig_info)

tadata.shape

df = tadata.obs[['fig_label','outlier']].value_counts().reset_index(name = 'count').pivot(index = 'fig_label', columns = 'outlier', values = 'count')
df2 = df.div(df.sum(axis = 1), axis = 0)*100

#### Making Figure 1 b,c)
cell_count = (adata.obs.fig_label.value_counts()).reset_index()
tadata_cell_count = (tadata.obs.fig_label.value_counts()).reset_index()

# Compute percentages
cell_count_perc = (adata.obs.fig_label.value_counts() / 12500 * 100).reset_index()
tadata_cell_count_perc = (tadata.obs.fig_label.value_counts() / 5000 * 100).reset_index()

# Rename columns for clarity
cell_count_perc.columns = ['fig_label', 'PARSE']
tadata_cell_count_perc.columns = ['fig_label', '10X']
cell_count.columns = ['fig_label', 'PARSE']
tadata_cell_count.columns = ['fig_label', '10X']

# Merge on 'fig_label' to align values correctly
cell_count_info = pd.merge(tadata_cell_count, cell_count, on='fig_label', how='outer')
cell_count_perc_info = pd.merge(tadata_cell_count_perc, cell_count_perc, on='fig_label', how='outer')

print(cell_count_info)
print(cell_count_perc_info)

import seaborn as sns
import matplotlib.pyplot as plt

fig, axes = plt.subplots(
    nrows=2, ncols=1,
    figsize=(4, 4),
    constrained_layout=True
)

# Melt the DataFrame to long format for easier plotting
cell_count_perc_melted = cell_count_perc_info.melt(id_vars='fig_label', var_name='Dataset', value_name='Percentage')
cell_count_melted = cell_count_info.melt(id_vars='fig_label', var_name='Dataset', value_name='Absolute')

# Create the bar plot
#plt.figure(figsize=(5, 4))
sns.barplot(data=cell_count_melted, x='fig_label', y='Absolute', hue='Dataset', palette = {'10X' :'navy', 'PARSE' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[0])
axes[0].set_xlabel('')
axes[0].set_ylabel('Cell count')
axes[0].legend(title="Dataset", bbox_to_anchor=(1.05, -0), loc='lower left')
axes[0].grid(False)
sns.barplot(data=cell_count_perc_melted, x='fig_label', y='Percentage', hue='Dataset', palette = {'10X' :'navy', 'PARSE' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[1])
axes[1].axhline(100, color = 'gold', linestyle = '-', linewidth = 2)
axes[1].set_xlabel("Patients")
axes[1].set_ylabel("Percentage of Cells (%)")
axes[1].legend_.remove()
axes[1].grid(False)
fig.align_ylabels()

plt.show()


#### Loading sequencing QC
reads = pd.read_csv('/Path/to/folder/with/QC_numbers.csv', sep = ';')

# Define target values
target_10X = 100000000
target_Parse = 250000000

# Function to calculate percentage based on Sample type
def calculate_percentage(row):
    if "10X" in row["Sample"]:  # For 10X samples
        return (row["Reads_all"] / target_10X) * 100
    elif "Parse" in row["Sample"]:  # For Parse samples
        return (row["Reads_all"] / target_Parse) * 100
    else:
        return None  # Handle unexpected cases

# Apply function to create the new column
reads_wt = reads[reads.Sample.str.contains('WT')]
reads_wt["percentage"] = reads_wt.apply(calculate_percentage, axis=1)
reads_wt['calc_mean_cell'] = reads_wt.Reads_assigned/reads_wt.Cell_count
reads_wt['add_mean_cell'] = reads_wt.Mean_read_per_cell-reads_wt.calc_mean_cell

# Display updated DataFrame
print(reads_wt)

#### Making figure 1 d,e)

# Count the number of 10X patients
num_10X = sum(patient.isdigit() for patient in reads_wt['Patient'].unique())
total_patients = len(reads_wt['Patient'].unique())

# Set xmax as the proportion of 10X patients in the total
xmax_10X = num_10X / total_patients

reads_wts = reads_wt[~reads_wt.Navn.str.contains('Sub')]

fig, axes = plt.subplots(
    nrows=2, ncols=2,
    figsize=(6, 4),
    constrained_layout=True
)

sns.barplot(data=reads_wts, x='Patient', y='Cell_count', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[0,0])
axes[0,0].axhline(12500, color = 'gold', linestyle = '-', linewidth = 2)
axes[0,0].axhline(5000, color = 'gold', linestyle = '-', linewidth = 2)
axes[0,0].set_xlabel('')
axes[0,0].set_ylabel('Cell count')
axes[0,0].legend_.remove()
#axes[0,0].text(- 0.2, 1.05,"b)", transform=axes[0,0].transAxes, fontsize=12, fontweight='bold', ha='right', va='bottom')

sns.barplot(data=reads_wts, x='Patient', y='cell_count_perc', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[1,0])
axes[1,0].axhline(100, color = 'gold', linestyle = '-', linewidth = 2)
axes[1,0].set_xlabel("")
axes[1,0].set_ylabel("Target cells (%)")
axes[1,0].legend_.remove()
#axes[1,0].text(- 0.2, 1.05,"c)", transform=axes[1,0].transAxes, fontsize=12, fontweight='bold', ha='right', va='bottom')

sns.barplot(data=reads_wts, x='Patient', y='Reads_assigned', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax=axes[0,1])
axes[0,1].legend_.remove()
axes[0,1].set_xlabel('')
axes[0,1].set_ylabel('Assigned reads')
#axes[0,1].text(- 0.2, 1.05,"d)", transform=axes[0,1].transAxes, fontsize=12, fontweight='bold', ha='right', va='bottom')

sns.barplot(data=reads_wts, x='Patient', y='calc_mean_cell', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[1,1])
axes[1,1].axhline(20000, color = 'gold', linestyle = '-', linewidth = 2)
axes[1,1].set_xlabel("")
axes[1,1].set_ylabel("Assigned reads/cell")
axes[1,1].legend_.remove()
#axes[1,1].text(- 0.2, 1.05,"e)", transform=axes[1,1].transAxes, fontsize=12, fontweight='bold', ha='right', va='bottom')

#fig.align_ylabels()

plt.show()
fig.savefig('figures/cellcount_assignedreads.svg', dpi = 300)

#### Making Supplementary Figure 1 a-d)
# Count the number of 10X patients
num_10X = sum(patient.isdigit() for patient in reads_wt['Patient'].unique())
total_patients = len(reads_wt['Patient'].unique())

# Set xmax as the proportion of 10X patients in the total
xmax_10X = num_10X / total_patients

reads_wts = reads_wt[~reads_wt.Navn.str.contains('Sub')]

max_width_mm = 190
max_width_inches = max_width_mm / 25.4  # Convert mm to inches
col_width = max_width_inches / 4
row_height = 1.5  # Keep similar aspect ratio

# Set up the figure and grid
fig, axes = plt.subplots(
    nrows=1, ncols=4,
    figsize=(max_width_inches, row_height),
    constrained_layout=True
)

sns.barplot(data=reads_wt, x='Patient', y='Reads_all', hue='Sample', order=sorted(reads_wt['Patient'].unique()), palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax=axes[0], width = 1.6)

n_bars = len(axes[0].patches)
n_categories = len(sorted(reads_wt['Patient'].unique()))
bars_per_category = n_bars // n_categories if n_categories > 0 else 0

# For each patient position
for i in range(n_categories):
    # Find the bars for this patient
    category_bars = [axes[0].patches[j] for j in range(n_bars) 
                    if j // bars_per_category == i]
    
    # For patients with only one bar, center it
    if len(category_bars) == 1:
        bar = category_bars[0]
        # Center the bar at integer position i
        bar.set_x(i - bar.get_width()/2)
axes[0].axhline(100000000, color = 'gold', linestyle = '-', linewidth = 2, xmax = xmax_10X)
axes[0].axhline(250000000, color = 'gold', linestyle = '-', linewidth = 2, xmin = xmax_10X)
axes[0].legend_.remove()
axes[0].set_ylabel('Total raw reads (1e8)')
axes[0].set_xlabel('')

sns.barplot(data=reads_wts, x='Patient', y='Valid_barcodes', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax=axes[1])
axes[1].legend_.remove()
axes[1].set_xlabel('')
axes[1].set_ylim(0,100)
axes[1].set_ylabel('Valid barcodes (%)')

sns.barplot(data=reads_wts, x='Patient', y='Transcriptome_VDJ_reads', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax=axes[2])
axes[2].legend_.remove()
axes[2].set_xlabel('')
axes[2].set_ylim(0,100)
axes[2].set_ylabel('Transcriptome reads (%)')

sns.barplot(data=reads_wts, x='Patient', y='sequencing_saturation', order=sorted(reads_wts['Patient'].unique()), hue='Sample', palette = {'10X-WT' :'navy', 'Parse-WT' :'firebrick'}, alpha = 0.8, edgecolor = 'black', ax=axes[3])
axes[3].legend_.remove()
axes[3].set_xlabel('')
axes[3].set_ylim(0,100)
axes[3].set_ylabel('Sequencing saturation (%)')


#fig.align_ylabels()
fig.savefig('figures/seqstats_supplementary.png', dpi = 300)
plt.show()

#### Making Figure 1 f,g)
dt_copy = dt.copy()
df_copy = df.copy()
df_copy['Source'] = '10X'
dt_copy['Source'] = 'PARSE'

combined_df = pd.concat([df_copy, dt_copy])

dt_copy = dt2.copy()
df_copy = df2.copy()
df_copy['Source'] = '10X'
dt_copy['Source'] = 'PARSE'

combined_df2 = pd.concat([df_copy, dt_copy])

fig, axes = plt.subplots(
    nrows=2, ncols=1,
    figsize=(3, 4),
    constrained_layout=True
)

sns.barplot(data=combined_df, x='fig_label', y=True, hue='Source', palette={'10X': 'navy', 'PARSE': 'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[0])
axes[0].set_xlabel("")
axes[0].set_ylabel("Outliers")
axes[0].legend_.remove()
axes[0].grid(False)
fig.align_ylabels()

sns.barplot(data=combined_df2, x='fig_label', y=True, hue='Source', palette={'10X': 'navy', 'PARSE': 'firebrick'}, alpha = 0.8, edgecolor = 'black', ax = axes[1])
axes[1].set_xlabel("Patients")
axes[1].set_ylabel("Outliers (%)")
axes[1].legend_.remove()
axes[1].grid(False)
fig.align_ylabels()
fig.savefig('figures/Outliers.svg', dpi = 300)

#Calculating Wilcoxon
# Merge dataframes on patient
merged_df = pd.concat([df, dt], keys = ['10X', 'PARSE'])

print(merged_df)

# Extract True counts and totals for each method
true_10X = merged_df.loc['10X', True]
true_PARSE = merged_df.loc['PARSE', True]

false_10X = merged_df.loc['10X', False]
false_PARSE = merged_df.loc['PARSE', False]

# Calculate proportions
prop_10X = true_10X / (true_10X + false_10X)
prop_PARSE = true_PARSE / (true_PARSE + false_PARSE)

# Run Wilcoxon on proportions
from scipy.stats import wilcoxon
statistic, p_value = wilcoxon(prop_10X, prop_PARSE)
print(f"Wilcoxon statistic: {statistic}")
print(f"p-value: {p_value}")


#Making Supplementary figure 1 e-l)
# Define custom y-axis labels for each variable
ylabels = {
    "log1p_total_counts": "Total Counts (log1p)",
    "log1p_n_genes_by_counts": "Number of Genes by Counts (log1p)",
    "pct_counts_in_top_20_genes": "Percentage Counts in Top 20 Genes",
    "pct_counts_mt": "Percentage Counts MT"
}

# Set font sizes globally
plt.rcParams.update({'font.size': 5})

# Variables and methods
variables = ['log1p_total_counts', "log1p_n_genes_by_counts", "pct_counts_in_top_20_genes", "pct_counts_mt"]
methods = ["10X", "PARSE"]  # Replace with actual data objects like `tadata` and `adata`

# Calculate figure size to not exceed 190mm width
max_width_mm = 190
max_width_inches = max_width_mm / 25.4  # Convert mm to inches
col_width = max_width_inches / len(variables)
row_height = 1.5  # Keep similar aspect ratio

# Set up the figure and grid
fig, axes = plt.subplots(
    nrows=len(methods), ncols=len(variables),
    figsize=(max_width_inches, row_height * len(methods)),
    constrained_layout=True
)

# Loop through methods and variables
for row, method in enumerate(methods):
    data = tadata if method == "10X" else adata  # Replace with actual data objects
    for col, var in enumerate(variables):
        ax = axes[row, col]
        # Violin plot
        sc.pl.violin(
            data,
            var,
            groupby="fig_label",
            ax=ax,
            show=False,
            size = 0.2
        )
        # Calculate thresholds
        if method == '10X':
            fig_labels = sorted(data.obs['fig_label'].unique())
            for idx, label in enumerate(fig_labels):
                group = data[data.obs.fig_label == label]
                if var == "pct_counts_mt":
                    ax.axhline(y=10, color='red', linestyle='--', linewidth=1, label = 'Threshold = 8')
                else:
                    low = np.median(group.obs[var]) - 5 * median_abs_deviation(group.obs[var])
                    up = np.median(group.obs[var]) + 5 * median_abs_deviation(group.obs[var])
                # Add hlines for thresholds
                    ax.hlines(y=up, xmin=idx - 0.25, xmax=idx + 0.25, color='black', linewidth=1)
                    ax.hlines(y=low, xmin=idx - 0.25, xmax=idx + 0.25, color='black', linewidth=1)
        else:
            if var == "pct_counts_mt":
                ax.axhline(y=10, color='red', linestyle='--', linewidth=1, label = 'Threshold = 8')
                #ax.legend()
            else:
                low = np.median(data.obs[var]) - 5 * median_abs_deviation(data.obs[var])
                up = np.median(data.obs[var]) + 5 * median_abs_deviation(data.obs[var])
            # Add hlines for thresholds (this assumes data doesn't have separate groups)
                ax.hlines(y=up, xmin=-1, xmax=7, color='black', linewidth=1)
                ax.hlines(y=low, xmin=-1, xmax=7, color='black', linewidth=1)
   
        # Set y-axis labels and remove titles
        if col == 0:
            ax.set_ylabel(f"{method}\n{ylabels[var]}", fontsize=5)
        else:
            ax.set_ylabel(ylabels[var], fontsize=5)
        
        # Remove title and set x-axis label
        ax.set_title('')
        if row == len(methods) - 1:  # Only bottom row gets x-axis labels
            ax.set_xlabel('Patients', fontsize=5)
        else:
            ax.set_xlabel('')
        
        # Set tick label font sizes
        ax.tick_params(axis='both', which='major', labelsize=5)

# Display combined plots
plt.grid(False)
plt.savefig('figures/QCparameters.png', dpi=300, bbox_inches='tight')
plt.show()


#Filtering the data based on outliers
pre_adata = adata

print(f"Total number of cells: {adata.n_obs}")
adata = adata[(~adata.obs.outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

pre_tadata = tadata

print(f"Total number of cells: {tadata.n_obs}")
tadata = tadata[(~tadata.obs.outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {tadata.n_obs}")


#Adding TCR sequence to Parse data
df_tcr_raw = pd.read_csv("/Path/to/folder/with/Parse/barcode_report.tsv", sep="\t")

df_tcr_raw.columns

cols_to_join = ['TRA_V', 'TRA_D', 'TRA_J', 'TRA_C', 'TRA_cdr3_aa',
       'TRA_read_count', 'TRA_transcript_count', 'TRB_V', 'TRB_D', 'TRB_J',
       'TRB_C', 'TRB_cdr3_aa', 'TRB_read_count', 'TRB_transcript_count',
       'secondary_TRA_V', 'secondary_TRA_D', 'secondary_TRA_J',
       'secondary_TRA_C', 'secondary_TRA_cdr3_aa', 'secondary_TRA_read_count',
       'secondary_TRA_transcript_count', 'secondary_TRB_V', 'secondary_TRB_D',
       'secondary_TRB_J', 'secondary_TRB_C', 'secondary_TRB_cdr3_aa',
       'secondary_TRB_read_count', 'secondary_TRB_transcript_count', 'isMultiplet']

df_tcr_raw.set_index('Barcode', inplace = True)
df_tcr_raw.index.name = None
adata.obs = adata.obs.join(df_tcr_raw[cols_to_join], how = "left")
adata.obs_names_make_unique()

adata.obs['has_TRB'] = adata.obs.TRB_cdr3_aa.notna()
adata.obs['has_TRA'] = adata.obs.TRA_cdr3_aa.notna()

#Filtering cells with no TCR
tcell_adata = adata[(adata.obs.has_TRB == True) | (adata.obs.has_TRA == True)]
print(tcell_adata.shape)
print(adata.shape)

tadata.obs['has_TCR'] = tadata.obs.cdr3s_aa.notna()
tcell_tadata = tadata[(tadata.obs.has_TCR == True)]
print(tcell_tadata.shape)
print(tadata.shape)

#Saving the filtered data
adata.write_h5ad('/Path/to/folder/PARSE_postQC.h5ad')
tcell_adata.write_h5ad('/Path/to/folder/PARSE_postQC_Tcells.h5ad')
tadata.write_h5ad('/Path/to/folder/10X_postQC.h5ad')
tcell_tadata.write_h5ad('/Path/to/folder/10X_postQC_Tcells.h5ad')
