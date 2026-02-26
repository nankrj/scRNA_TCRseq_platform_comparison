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
from pyensembl import EnsemblRelease
import re
from scipy.stats import pearsonr
from scipy.stats import gaussian_kde
from collections import defaultdict
from sklearn.linear_model import LinearRegression

ens98 = EnsemblRelease(release=98, species='human')
ens98.download()
ens98.index()

import json

# Load both
with open('/faststorage/project/TCRseq_AZ/Analysis/Nanna/SingleCell/Transcript_length/transcript_lengths_ensembl98.json', 'r') as f:
    lengths_98 = json.load(f)

with open('/faststorage/project/TCRseq_AZ/Analysis/Nanna/SingleCell/Transcript_length/transcript_lengths_ensembl109.json', 'r') as f:
    lengths_109 = json.load(f)

genes_in_both_ensembls = (set(lengths_98.keys()) & set(lengths_109.keys()))

# Load data
tadata_all = sc.read_h5ad('/Path/to/folder/with/10X_adata_postQC.h5ad')
adata_all = sc.read_h5ad('/Path/to/folder/with/PARSE_adata_postQC.h5ad')

def get_ensembl_id(gene):
    try:
        genes = ens98.genes_by_name(gene)  # Get list of gene objects
        if genes:
            return genes[0].gene_id  # Return the first Ensembl ID
    except:
        genes = ens98.genes_by_name(re.search(r'^(.*?)-1', gene).group(1))
        if genes:
            return genes[1].gene_id  # Return the first Ensembl ID
        #return "Unknown"  # If not found, return "Unknown"

# Apply function to get Ensembl IDs
tadata_all.var["Ensemble_ID"] = [get_ensembl_id(gene) for gene in tadata_all.var_names]

# Create figure with 4 subplots
fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)

# Extract gene names
P_genes_all = set(adata_all.var.gene_id) & genes_in_both_ensembls
T_genes_all = set(tadata_all.var.Ensemble_ID) & genes_in_both_ensembls

# Function to create Venn diagram with absolute counts + percentages
def create_venn(ax, set1, set2, title):
    venn = venn2(subsets=(1, 1, 1), set_labels=('10X', 'PARSE'), set_colors=('navy', 'firebrick'), alpha=0.72, ax=ax)

    # Calculate values
    only_set1 = len(set1 - set2)
    only_set2 = len(set2 - set1)
    intersection = len(set1 & set2)
    total = only_set1 + only_set2 + intersection

    # Format labels with absolute counts and percentages
    labels = {
        '10': f"{only_set1}\n({only_set1/len(set1):.1%})",
        '01': f"{only_set2}\n({only_set2/len(set2):.1%})",
        '11': f"{intersection}"
    }

    # Update text labels
    for key, text in labels.items():
        label = venn.get_label_by_id(key)
        if label:
            label.set_text(text)
            label.set_fontsize(12)  # Adjust font size

    for patch in venn.patches:
        if patch:
            patch.set_edgecolor("black")  # Set edge color
            patch.set_linewidth(1)  # Adjust line width

    ax.set_title(title, pad=8)

# Venn diagram for all genes
create_venn(axes[0], T_genes_all, P_genes_all, 'All genes')

# Filter genes present in >20 cells
adata = adata_all.copy()
sc.pp.filter_genes(adata, min_cells=20)
tadata = tadata_all.copy()
sc.pp.filter_genes(tadata, min_cells=20)

P_genes = set(adata.var.gene_id) & genes_in_both_ensembls
T_genes = set(tadata.var.Ensemble_ID) & genes_in_both_ensembls

# Venn diagram for genes present in >20 cells
create_venn(axes[1], T_genes, P_genes, 'Genes present in >20 cells')

# Normalize and find highly variable genes
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)

sc.pp.normalize_total(tadata, target_sum=1e4)
sc.pp.log1p(tadata)
sc.pp.highly_variable_genes(tadata, min_mean=0.0125, max_mean=3, min_disp=0.25)

hv_genes_P = set(adata.var[adata.var.highly_variable].gene_id) & genes_in_both_ensembls
hv_genes_X = set(tadata.var[tadata.var.highly_variable].Ensemble_ID) & genes_in_both_ensembls

# Venn diagram for highly variable genes
create_venn(axes[2], hv_genes_X, hv_genes_P, 'Highly variable genes')

# Find top 1000 genes by total counts
sc.pp.calculate_qc_metrics(tadata, percent_top=[20], log1p=True, inplace=True)

P1000_genes = set(adata.var.loc[adata.var.total_counts.nlargest(1005).index, "gene_id"]) & genes_in_both_ensembls
T1000_genes = set(tadata.var.loc[tadata.var.total_counts.nlargest(1002).index, "Ensemble_ID"]) & genes_in_both_ensembls

# Venn diagram for top 1000 genes
create_venn(axes[3], T1000_genes, P1000_genes, 'Top 1000 genes')

plt.savefig('Parse_10X/figures/common_genes_revision.svg', dpi = 300)

ens109 = EnsemblRelease(release=109, species='human')
ens109.download()
ens109.index()

common_1000 = P1000_genes & T1000_genes
T_unique_1000 = T1000_genes - P1000_genes
P_unique_1000 = P1000_genes - T1000_genes

common_hv = hv_genes_P & hv_genes_X
T_unique_hv = hv_genes_X - hv_genes_P
P_unique_hv = hv_genes_P - hv_genes_X

common = P_genes & T_genes
T_unique = T_genes - P_genes
P_unique = P_genes - T_genes

common_all = P_genes_all & T_genes_all
T_unique_all = T_genes_all - P_genes_all
P_unique_all = P_genes_all - T_genes_all

all_genes = common_all | T_unique_all | P_unique_all

print(f'Percent of genes remaining in the 10X data after filtering for >20 : {len(T_genes)/len(T_genes_all)}')
print(f'Percent of genes remaining in the Parse data after filtering for >20 : {len(P_genes)/len(P_genes_all)}')
print(f'Percent of Parse unique genes remaining after filtering for >20 : {len(P_unique_all & P_unique)/len(P_unique_all)}')
print(f'Percent of 10X or Parse unique genes in highly variable or top category that are detected by both methods originally : {len((P_unique_hv | T_unique_hv | P_unique_1000 | T_unique_1000) & common_all)/len(P_unique_hv | T_unique_hv | P_unique_1000 | T_unique_1000)}')

import matplotlib.pyplot as plt
import json
from scipy import stats
import numpy as np

# Load the transcript length dictionaries
with open('/faststorage/project/TCRseq_AZ/Analysis/Nanna/SingleCell/Transcript_length/transcript_lengths_ensembl109.json', 'r') as f:
    lengths_109 = json.load(f)

# Create subsets from gene sets
def get_lengths(gene_set, ensembl_dict):
    """Get transcript lengths for a set of genes"""
    return {gene: ensembl_dict[gene] for gene in gene_set if gene in ensembl_dict}

datasets = [
    {
        'name': 'All transcripts',
        'tenx': get_lengths(T_unique_all, lengths_109),
        'parse': get_lengths(P_unique_all, lengths_109),
        'common': get_lengths(common_all, lengths_109),
    },
    {
        'name': '>20 cells',
        'tenx': get_lengths(T_unique, lengths_109),
        'parse': get_lengths(P_unique, lengths_109),
        'common': get_lengths(common, lengths_109),
    },
    {
        'name': 'Highly variable',
        'tenx': get_lengths(T_unique_hv, lengths_109),
        'parse': get_lengths(P_unique_hv, lengths_109),
        'common': get_lengths(common_hv, lengths_109),
    },
    {
        'name': 'Top 1000',
        'tenx': get_lengths(T_unique_1000, lengths_109),
        'parse': get_lengths(P_unique_1000, lengths_109),
        'common': get_lengths(common_1000, lengths_109),
    }
]

# Create boxplot
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
colors = ['navy', 'green', 'firebrick']

def add_significance_bracket(ax, x1, x2, y, p_val, median1, median2):
    """Add a significance bracket with fold-change"""
    # Skip if either median is None or zero
    if median1 is None or median2 is None or median1 == 0 or median2 == 0:
        return
        
    fold_change = median2 / median1 if median1 > 0 else 0
    
    # Format label
    if fold_change > 1:
        label = f'{fold_change:.1f}x'
    else:
        label = f'{1/fold_change:.1f}x'
    
    # Draw the bracket
    ax.plot([x1, x1, x2, x2], [y, y*1.05, y*1.05, y], 'k-', lw=1.5)
    # Add the fold-change label
    ax.text((x1 + x2) / 2, y*1.15, label, ha='center', fontsize=9)

for idx, dataset in enumerate(datasets):
    ax = axes[idx]
    
    tenx_vals = list(dataset['tenx'].values())
    parse_vals = list(dataset['parse'].values())
    common_vals = list(dataset['common'].values())
    
    data_to_plot = [tenx_vals, common_vals, parse_vals]
    labels = [f'10X-unique\n(n={len(tenx_vals)})', 
              f'Common\n(n={len(common_vals)})',
              f'Parse-unique\n(n={len(parse_vals)})']
    
    bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
    
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    
    ax.set_ylabel('Transcript Length (bp - log10)')
    ax.set_yscale('log')
    ax.set_title(dataset['name'])
    ax.grid(axis='y', alpha=0.3, which='both')
    
    # Pairwise Mann-Whitney U tests
    # Pairwise Mann-Whitney U tests (only if both groups have data)
    p_10x_parse = np.nan
    p_10x_common = np.nan
    p_parse_common = np.nan
    
    if len(tenx_vals) > 0 and len(parse_vals) > 0:
        u_stat_10x_parse, p_10x_parse = stats.mannwhitneyu(tenx_vals, parse_vals)
    
    if len(tenx_vals) > 0 and len(common_vals) > 0:
        u_stat_10x_common, p_10x_common = stats.mannwhitneyu(tenx_vals, common_vals)
    
    if len(common_vals) > 0 and len(parse_vals) > 0:
        u_stat_parse_common, p_parse_common = stats.mannwhitneyu(parse_vals, common_vals)
        
    # Get y-axis range for bracket placement
    y_min, y_max = ax.get_ylim()
    y_range = np.log10(y_max) - np.log10(y_min)
    bracket_y = 10 ** (np.log10(y_max) - y_range * 0.04)
    
    # Add significance brackets only if all groups have data
    if tenx_vals and common_vals:
        add_significance_bracket(ax, 1, 2, bracket_y, p_10x_common, np.median(tenx_vals), np.median(common_vals))
    
    if common_vals and parse_vals:
        add_significance_bracket(ax, 2, 3, bracket_y * 1.2, p_parse_common, np.median(common_vals), np.median(parse_vals))
    
    if tenx_vals and parse_vals:
        add_significance_bracket(ax, 1, 3, bracket_y * 2.4, p_10x_parse, np.median(tenx_vals), np.median(parse_vals))
    
    # Print stats for this dataset
    print(f"\n=== {dataset['name']} ===")
    print(f"10X-unique: median={np.median(tenx_vals):.0f} (n={len(tenx_vals)})")
    print(f"Parse-unique: median={np.median(parse_vals):.0f} (n={len(parse_vals)})")
    print(f"Shared: median={np.median(common_vals):.0f} (n={len(common_vals)})")
    print(f"Mann-Whitney U tests:")
    print(f"  10X vs Parse: p={p_10x_parse:.2e}")
    print(f"  10X vs Shared: p={p_10x_common:.2e}")
    print(f"  Parse vs Shared: p={p_parse_common:.2e}")

plt.tight_layout()
plt.savefig('transcript_length.svg', dpi=300, bbox_inches='tight')
plt.show()

biotype_dict = {}
for i, gene in enumerate(all_genes):
    try:
        biotype_dict[gene] = ens109.gene_by_id(gene).biotype
    except:
        pass
    
    if (i + 1) % 1000 == 0:
        print(f"  Fetched {i + 1}/{len(all_genes)} genes")

print(f"Total genes with biotypes: {len(biotype_dict)}\n")

def categorize_biotype(biotype):
    """Standardize biotype names"""
    biotype_lower = biotype.lower()
    
    if 'pseudogene' in biotype_lower:
        return 'pseudogene'
    elif 'ig_' in biotype_lower:
        return 'IG_x_gene'
    elif 'tr_' in biotype_lower:
        return 'TR_x_gene'
    elif biotype in ['miRNA', 'snRNA', 'snoRNA', 'misc_RNA', 'sRNA', 'scaRNA', 'rRNA']:
        return 'small non-coding RNAs'
    elif biotype in ['vault_RNA', 'scRNA', 'ribozyme', 'Mt_rRNA', 'Mt_tRNA', 'artifact']:
        return 'Other RNAs'
    else:
        return biotype

def count_biotypes(gene_set, biotype_dict):
    """
    Count biotypes in a gene set using pre-fetched biotype dictionary
    """
    counts = defaultdict(int)
    for gene in gene_set:
        if gene in biotype_dict:
            biotype = categorize_biotype(biotype_dict[gene])
            counts[biotype] += 1
    return counts

results = {
    'Common; all': count_biotypes(common_all, biotype_dict),
    'Common; >20 cells': count_biotypes(common, biotype_dict),
    'Common; Highly variable': count_biotypes(common_hv, biotype_dict),
    'Common; top 1000': count_biotypes(common_1000, biotype_dict),
    '10X; all': count_biotypes(T_unique_all, biotype_dict),
    '10X; >20 cells': count_biotypes(T_unique, biotype_dict),
    '10X; Highly variable': count_biotypes(T_unique_hv, biotype_dict),
    '10X; top 1000': count_biotypes(T_unique_1000, biotype_dict),
    'Parse; all': count_biotypes(P_unique_all, biotype_dict),
    'Parse; >20 cells': count_biotypes(P_unique, biotype_dict),
    'Parse; Highly variable': count_biotypes(P_unique_hv, biotype_dict),
    'Parse; top 1000': count_biotypes(P_unique_1000, biotype_dict),
}

combined_df = pd.DataFrame(results).fillna(0).astype(int)

biotype_order = [
    'protein_coding',
    'lncRNA',
    'small non-coding RNAs',
    'pseudogene',
    'IG_x_gene',
    'TR_x_gene', 
    'TEC',
    'Other RNAs'
]
biotype_order = [b for b in biotype_order if b in combined_df.index]
combined_df = combined_df.loc[biotype_order]

print(combined_df)

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

sns.set_theme(style="white")

def plot_stacked_barplot_subplot(ax, df, color_dict, legend=False):
    """Create a stacked barplot on a given axis"""
    biotypes = df.index
    categories = df.columns
    colors = [color_dict.get(biotype, "#999999") for biotype in biotypes]
    
    bottoms = np.zeros(len(categories))
    
    for i, biotype in enumerate(biotypes):
        values = df.loc[biotype].values
        ax.bar(categories, values, bottom=bottoms, width=0.85,
               label=biotype, color=colors[i], edgecolor='black', alpha=0.8)
        bottoms += values
    
    # Extract clean labels from column names
    clean_labels = []
    for col in categories:
        if "10X" in col:
            clean_labels.append("10X")
        elif "Common" in col:
            clean_labels.append("Common")
        elif "Parse" in col or "PARSE" in col:
            clean_labels.append("Parse")
        else:
            clean_labels.append(col)
    
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(clean_labels, rotation=45, ha='right')
    ax.set_ylabel("Percentage (%)")
    
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, title="Biotype", bbox_to_anchor=(1.05, 1), loc="upper left")

# Define biotype order and properties
biotype_order = [
    'protein_coding',
    'lncRNA',
    'small non-coding RNAs',
    'pseudogene',
    'IG_x_gene',
    'TR_x_gene', 
    'TEC',
    'Other RNAs'
]

biotype_colors = {
    'protein_coding': 'tab:blue',
    'lncRNA': 'tab:orange',
    'small non-coding RNAs': 'tab:green',
    'pseudogene': 'tab:red',
    'IG_x_gene': 'tab:purple',
    'TR_x_gene': 'tab:brown', 
    'TEC': 'tab:pink',
    'Other RNAs': 'tab:gray'
}

# Prepare data
combined_percent = combined_df.div(combined_df.sum(axis=0), axis=1) * 100

# Define column order: 10X, Common, Parse (for each category)
categories_main = ['all', '>20 cells', 'Highly variable', 'top 1000']
platforms = ['10X', 'Common', 'Parse']

# Create 1x4 subplots
fig, axes = plt.subplots(1, 4, figsize=(18, 4))

for subplot_idx, category in enumerate(categories_main):
    ax = axes[subplot_idx]
    
    # Get columns for this category
    category_columns = [col for col in combined_percent.columns 
                       if category in col]
    category_data = combined_percent[category_columns]
    
    # Reorder by platform (10X, Common, Parse)
    ordered_cols = []
    for platform in platforms:
        matching = [col for col in category_columns if platform in col]
        ordered_cols.extend(matching)
    
    category_data = category_data[ordered_cols]
    
    # Ensure only existing biotypes are included
    biotype_order_filtered = [b for b in biotype_order if b in category_data.index]
    category_data = category_data.loc[biotype_order_filtered]
    
    # Plot on this subplot
    plot_stacked_barplot_subplot(ax, category_data, biotype_colors, legend=(subplot_idx == 3))
    ax.set_title(category, fontweight='bold')

plt.tight_layout()
plt.savefig('biotype_distribution.svg', dpi=300, bbox_inches='tight')
plt.show()

fig, axes = plt.subplots(
    nrows=1, ncols=4,
    figsize=(14, 3),
    constrained_layout=True
)
main_categories = ["All genes", "Genes present in >20 cells", "Highly variable genes", "Top 1000 genes"]
genes_set1 = [P_genes_all, P_genes, hv_genes_P, P1000_genes]  # 10X genes
genes_set2 = [T_genes_all, T_genes, hv_genes_X, T1000_genes]  # PARSE genes
genes = [
    (genes_set1[0] | genes_set2[0]),
    (genes_set1[1] | genes_set2[1]),
    (genes_set1[2] | genes_set2[2]),
    (genes_set1[3] | genes_set2[3])
]
# Define your custom colors
custom_colors = ['firebrick', 'navy', 'green']  # 10X only, PARSE only, both '#a15898'
# Create a custom colormap
from matplotlib.colors import ListedColormap
custom_cmap = ListedColormap(custom_colors)
# Flatten the axes array to make indexing easier
axes = axes.flatten()
for idx, g in enumerate(genes):
    d1 = adata_all[:, adata_all.var.gene_id.isin(g)]
    d2 = tadata_all[:, tadata_all.var.Ensemble_ID.isin(g)]
    
    # Check if there are samples in both datasets for this label
    if d1.shape[0] > 0 and d2.shape[0] > 0:
        parse_avg = pd.DataFrame(d1.X.mean(axis=0).A1, index=d1.var.gene_id, columns=["Parse"])
        tenx_avg = pd.DataFrame(d2.X.mean(axis=0).A1, index=d2.var.Ensemble_ID, columns=["10X"])
        merged = parse_avg.join(tenx_avg, how="outer").fillna(0)
        
        # Create a color array based on gene membership
        color_array = []
        for gene in merged.index:
            if gene in genes_set1[idx] and gene in genes_set2[idx]:
                color_array.append(2)  # Both sets (purple: #a15898)
            elif gene in genes_set1[idx]:
                color_array.append(0)  # 10X only (navy)
            elif gene in genes_set2[idx]:
                color_array.append(1)  # PARSE only (firebrick)
            else:
                color_array.append(0)  # Default (should not happen)
        
        # Calculate Pearson correlation for non-zero values
        # We'll use the log-transformed values for correlation to match the plot
        log_10x = np.log10(merged["10X"] + 1e-6)
        log_parse = np.log10(merged["Parse"] + 1e-6)
        
        # Calculate correlation
        corr, p_value = pearsonr(log_10x, log_parse)
        
        # Plot with pseudo-log scale and color by set membership
        sc = axes[idx].scatter(
            log_10x,
            log_parse,
            c=color_array,
            cmap=custom_cmap,  # Using your custom colormap
            alpha=0.3,  # Make points more transparent
            s=1,  # Make points smaller to avoid overcrowding
            rasterized=True  # For better performance with many points
        )
        
        # Add grid lines to help with visual reference
        axes[idx].grid(True, linestyle='--', alpha=0.3)
        
        # Add a diagonal line to show y=x
        max_val = max(np.log10(merged["10X"].max() + 1e-6), np.log10(merged["Parse"].max() + 1e-6))
        min_val = min(np.log10(merged["10X"].min() + 1e-6), np.log10(merged["Parse"].min() + 1e-6))
        axes[idx].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)

        axes[idx].set_title(f"{main_categories[idx]}")
        axes[idx].set_xlabel('10X Average Gene Expression (log10)')
        axes[idx].set_ylim(top = 2.7)
        axes[idx].set_xlim(right = 2.7)
        
               
        # Add a legend to the first subplot and include correlation coefficient and intercept
        if idx == 0:
            axes[idx].set_ylabel('Parse Average Gene Expression (log10)')
            # Add linear regression line
            X = log_10x.values.reshape(-1, 1)
            y = log_parse.values
            
            # Create and fit the linear regression model
            model = LinearRegression()
            model.fit(X, y)
            
            # Get the slope and intercept
            slope = model.coef_[0]
            log_intercept = model.intercept_
        
            # Calculate the original scale intercept (10^intercept)
            # This gives the Parse value when 10X = 1 (since 10^0 = 1)
            original_intercept = 10 ** log_intercept - 1e-6  # Subtracting the pseudo-count we added
            
            # Create a range of x values for the regression line
            x_range = np.linspace(min_val, max_val, 100).reshape(-1, 1)
            y_pred = model.predict(x_range)
            
            # Plot the regression line in yellow
            axes[idx].plot(x_range, y_pred, color='black', linewidth=2)
            
            # Add text for correlation coefficient and intercept
            axes[idx].text(-5.7, 1.3, f'log10 intercept = {log_intercept:.3f}\noriginal scale = {original_intercept:.3f}')
            axes[idx].text(0,-5, f'r = {corr:.3f}')
        else:
            axes[idx].set_ylabel('')
            # Add intercept text to all other plots
            #axes[idx].text(-5, 2, f'intercept = {intercept:.3f}')
    else:
        axes[idx].text(0.5, 0.5, f"No data for {main_categories[idx]}", 
                      horizontalalignment='center', verticalalignment='center')
        axes[idx].set_title(main_categories[idx])
    
#plt.tight_layout()
plt.savefig('Avg_gene_exp_corr_revision.svg', bbox_inches='tight', dpi=300)

fig, axes = plt.subplots(
    nrows=2, ncols=4,
    figsize=(14, 6),
    constrained_layout=True
)

# Flatten the axes array to make indexing easier
axes = axes.flatten()



for idx, i in enumerate(['1', '2', '3', '4', '5', '6', '7']):
    d1 = adata_all[adata_all.obs.fig_label == i]
    d1x = d1[:, d1.var.gene_id.isin(all_genes)]
    d2 = tadata_all[tadata_all.obs.fig_label == i]
    d2x = d2[:, d2.var.Ensemble_ID.isin(all_genes)]
    
    # Check if there are samples in both datasets for this label
    if d1x.shape[0] > 0 and d2x.shape[0] > 0:
        parse_avg = pd.DataFrame(d1x.X.mean(axis=0).A1, index=d1x.var.gene_id, columns=["Parse"])
        tenx_avg = pd.DataFrame(d2x.X.mean(axis=0).A1, index=d2x.var.Ensemble_ID, columns=["10X"])
        merged = parse_avg.join(tenx_avg, how="outer").fillna(0)
        
        # Add a small constant to avoid log(0) issues when plotting
        #epsilon = 1e-6  # Smaller than any meaningful expression value
        log_10x = np.log10(merged["10X"] + 1e-6)
        log_parse = np.log10(merged["Parse"] + 1e-6)

        # Calculate correlation
        corr, p_value = pearsonr(log_10x, log_parse)
        
        # Plot with pseudo-log scale (log1p is natural log of 1+x)
        sc = axes[idx].scatter(
            log_10x,
            log_parse,
            alpha=0.3,  # Make points more transparent
            s=1,  # Make points smaller to avoid overcrowding
            rasterized=True  # For better performance with many points
        )
        
        # Add grid lines to help with visual reference
        axes[idx].grid(True, linestyle='--', alpha=0.3)
        
        # Add a diagonal line to show y=x
        # Add a diagonal line to show y=x
        max_val = max(np.log10(merged["10X"].max() + 1e-6), np.log10(merged["Parse"].max() + 1e-6))
        min_val = min(np.log10(merged["10X"].min() + 1e-6), np.log10(merged["Parse"].min() + 1e-6))
        axes[idx].plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)

        # Add linear regression line
        X = log_10x.values.reshape(-1, 1)
        y = log_parse.values
        
        # Create and fit the linear regression model
        model = LinearRegression()
        model.fit(X, y)
        
        # Get the slope and intercept
        slope = model.coef_[0]
        log_intercept = model.intercept_
    
        # Calculate the original scale intercept (10^intercept)
        # This gives the Parse value when 10X = 1 (since 10^0 = 1)
        original_intercept = 10 ** log_intercept - 1e-6  # Subtracting the pseudo-count we added
        
        # Create a range of x values for the regression line
        x_range = np.linspace(min_val, max_val, 100).reshape(-1, 1)
        y_pred = model.predict(x_range)
        
        # Plot the regression line in yellow
        axes[idx].plot(x_range, y_pred, color='black', linewidth=2)
        
        # Add text for correlation coefficient and intercept
        axes[idx].text(-5.85, 1.3, f'log10 intercept = {log_intercept:.3f}\noriginal scale = {original_intercept:.3f}')
        axes[idx].text(0,-5.5, f'r = {corr:.3f}')
        
        axes[idx].set_title(f'Average Gene Expression (log scale)')
        axes[idx].set_xlabel(f'10X Patient {i}')
        axes[idx].set_ylabel(f'Parse Patient {i}')
        
    else:
        axes[idx].text(0.5, 0.5, f"No data for label {i}", 
                      horizontalalignment='center', verticalalignment='center')
        axes[idx].set_title(f'Fig Label {i}')
    

# Hide the unused subplot
axes[7].set_visible(False)

#plt.tight_layout()
plt.savefig('Avg_gene_exp_corr_pr_patient_revision.png', bbox_inches='tight', dpi=300)


# Function to create Bland-Altman plot
def bland_altman_plot(x, y, title=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 3))  # Fixed: was "fix = g, ax"
    
    # Convert to arrays for calculations
    x = np.array(x)
    y = np.array(y)
    
    # Calculate mean and difference between methods
    mean = (x + y) / 2
    diff = y - x  # Parse - 10X
    
    # Calculate mean difference and limits of agreement
    md = np.mean(diff)
    sd = np.std(diff, axis=0)
    upper_loa = md + 1.96 * sd
    lower_loa = md - 1.96 * sd
    
    # Create the plot
    ax.scatter(mean, diff, alpha=0.3, s=2)
    ax.axhline(md, color='k', linestyle='-', linewidth=1)
    ax.axhline(upper_loa, color='r', linestyle='--', linewidth=1)
    ax.axhline(lower_loa, color='r', linestyle='--', linewidth=1)
    
    # Add text annotations
    ax.text(np.max(mean) + 0.5, md, f'Mean diff: {md:.3f}', 
            verticalalignment='center', fontsize=10)
    ax.text(np.max(mean) + 0.5, upper_loa, f'+1.96 SD: {upper_loa:.3f}', 
            verticalalignment='center', fontsize=10)
    ax.text(np.max(mean) + 0.5, lower_loa, f'-1.96 SD: {lower_loa:.3f}', 
            verticalalignment='center', fontsize=10)
    
    # Set labels and title
    ax.set_xlabel('Mean of 10X and Parse measurements')
    ax.set_ylabel('Difference (Parse - 10X)')
    ax.set_title(title)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    return fig, ax  # Return both figure and axis

# Prepare data
parse_avg = pd.DataFrame(d1x.X.mean(axis=0).A1, index=d1x.var.gene_id, columns=["Parse"])
tenx_avg = pd.DataFrame(d2x.X.mean(axis=0).A1, index=d2x.var.Ensemble_ID, columns=["10X"])
merged = parse_avg.join(tenx_avg, how="outer").fillna(0)

# Use log-transformed values
log_10x = np.log10(merged["10X"] + 1e-6)
log_parse = np.log10(merged["Parse"] + 1e-6)

# Create Bland-Altman plot
fig, ax = bland_altman_plot(log_10x, log_parse, title='Gene Expression Comparison')

# Save the figure
fig.savefig('Bland_Altman_gene_exp_allgenes_revision.png', bbox_inches='tight', dpi=300)

plt.show()

