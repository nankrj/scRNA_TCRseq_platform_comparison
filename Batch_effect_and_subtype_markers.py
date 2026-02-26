import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import json
import scipy.sparse

adata = sc.read_h5ad('/Path/to//PARSE_postQC_Tcells.h5ad')
latent = np.load('/Path/to//latent_space.npy',allow_pickle=True)
adata.obsm['X_scVI'] = latent

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pp.neighbors(adata, use_rep='X_pca', n_pcs=30)
sc.tl.umap(adata)

tadata = sc.read_h5ad('/Path/to/10X_postQC_Tcells.h5ad')
latent = np.load('/Path/to/latent_space.npy',allow_pickle=True)
tadata.obsm['X_scVI'] = latent

sc.tl.pca(tadata)
sc.pl.pca_variance_ratio(tadata, log=True, n_pcs=50)
sc.pp.neighbors(tadata, use_rep='X_pca', n_pcs=30)
sc.tl.umap(tadata)

#### Figures 3a-d

fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(2, 4))

# Plot 10X data on top
sc.pl.umap(tadata, color='patient', size=2, ax=ax_top, show=False, title='10X', legend_loc = None)

# Plot Parse data on bottom  
sc.pl.umap(adata, color='patient', size=2, ax=ax_bottom, show=False, title='Parse', legend_loc = None)

plt.tight_layout()
plt.savefig('prebatch_patient_umap.png')


sc.pp.neighbors(adata, use_rep='X_scVI')
sc.tl.umap(adata)

sc.pp.neighbors(tadata, use_rep='X_scVI')
sc.tl.umap(tadata)

# Plot 10X data on top
sc.pl.umap(tadata, color='patient', size=2, ax=ax_top, show=False, title='10X', legend_loc = None)

# Plot Parse data on bottom  
sc.pl.umap(adata, color='patient', size=2, ax=ax_bottom, show=False, title='Parse', legend_loc = None)

plt.tight_layout()
plt.savefig('postbatch_patient_umap.png')

#### Supplementary Figure 2a

#Data retrieved from batch correction runfile output
#From 10X
x = [5.96e+3,5.96e+3,5.96e+3,4.86e+3,4.86e+3,4.86e+3,4.8e+3,4.8e+3,4.8e+3,4.77e+3,4.77e+3,4.77e+3,4.75e+3,4.75e+3,4.75e+3,4.73e+3,4.73e+3,4.73e+3,4.71e+3,4.71e+3,4.71e+3,4.7e+3,4.7e+3,4.7e+3,4.68e+3,4.68e+3,4.68e+3,
     4.67e+3,4.67e+3,4.67e+3,4.66e+3,4.66e+3,4.66e+3,4.65e+3,4.65e+3,4.65e+3,4.64e+3,4.64e+3,4.64e+3,4.63e+3,4.63e+3,4.63e+3,4.62e+3,4.62e+3,4.62e+3,4.62e+3,4.62e+3,4.62e+3,4.61e+3,4.61e+3,4.61e+3,4.6e+3,4.6e+3,4.6e+3,
     4.59e+3,4.59e+3,4.59e+3,4.59e+3,4.59e+3,4.59e+3,4.58e+3,4.58e+3,4.58e+3,4.58e+3,4.58e+3,4.58e+3,4.57e+3,4.57e+3,4.57e+3,4.57e+3,4.57e+3,4.57e+3,4.56e+3,4.56e+3,4.56e+3,4.56e+3,4.56e+3,4.56e+3,4.55e+3,4.55e+3,
     4.55e+3,4.55e+3,4.55e+3,4.55e+3,4.54e+3,4.54e+3,4.54e+3,4.54e+3,4.54e+3,4.54e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.53e+3,4.52e+3,4.52e+3,4.52e+3,4.52e+3,4.52e+3,4.52e+3,4.52e+3,4.52e+3,
     4.52e+3,4.52e+3,4.52e+3,4.52e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.51e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.5e+3,4.49e+3,
     4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.49e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,
     4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.48e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,
     4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.47e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,
     4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,4.46e+3,
     4.46e+3,4.46e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,
     4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.45e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,
     4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,
     4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3,4.44e+3]

n = ['1','2','2','2','3','3','3','4','4','4','5','5','5','6','6','6','7','7','7','8','8','8','9','9','9','10','10','10','11','11','11','12','12','12','13','13','13','14','14','14','15','15','15','16','16','16','17',
     '17','17','18','18','18','19','19','19','20','20','20','21','21','21','22','22','22','23','23','23','24','24','24','25','25','25','26','26','26','27','27','27','28','28','28','29','29','29','30','30','30','31',
     '31','31','32','32','32','33','33','33','34','34','34','35','35','35','36','36','36','37','37','37','38','38','38','39','39','39','40','40','40','41','41','41','42','42','42','43','43','43','44','44','44','45',
     '45','45','46','46','46','47','47','47','48','48','48','49','49','49','50','50','50','51','51','51','52','52','52','53','53','53','54','54','54','55','55','55','56','56','56','57','57','57','58','58','58','59',
     '59','59','60','60','60','61','61','61','62','62','62','63','63','63','64','64','64','65','65','65','66','66','66','67','67','67','68','68','68','69','69','69','70','70','70','71','71','71','72','72','72','73',
     '73','73','74','74','74','75','75','75','76','76','76','77','77','77','78','78','78','79','79','79','80','80','80','81','81','81','82','82','82','83','83','83','84','84','84','85','85','85','86','86','86','87',
     '87','87','88','88','88','89','89','89','90','90','90','91','91','91','92','92','92','93','93','93','94','94','94','95','95','95','96','96','96','97','97','97','98','98','98','99','99','99','100','100','100','101',
     '101','101','102','102','102','103','103','103','104','104','104','105','105','105','106','106','106','107','107','107','108','108','108','109','109','109','110','110','110','111','111','111','111']

l = pd.DataFrame({'epoch': n, 'loss': x})
l.drop_duplicates()

#From Parse
parse_training = pd.DataFrame(
    {'epoch' : ['1','2','2','2','3','3','3','4','4','4','5','5','5','6','6','6','7','7','7','8','8','8','9','9','9','10','10','10','11','11','11','12','12','12','13','13','13','14','14','14','15','15','15',
                '16','16','16','17','17','17','18','18','18','19','19','19','20','20','20','21','21','21','22','22','22','23','23','23','24','24','24','25','25','25','26','26','26','27','27','27','28','28',
                '28','29','29','29','30','30','30','31','31','31','32','32','32','33','33','33','34','34','34','35','35','35','36','36','36','37','37','37','38','38','38','39','39','39','40','40','40','41',
                '41','41','42','42','42','43','43','43','44','44','44','45','45','45','46','46','46','47','47','47','48','48','48','49','49','49','50','50','50','51','51','51','52','52','52','53','53','53',
                '54','54','54','55','55','55','56','56','56','57','57','57','58','58','58','59','59','59','60','60','60','61','61','61','62','62','62','63','63','63','64','64','64','65','65','65','66','66',
                '66','67','67','67','68','68','68','69','69','69','70','70','70','71','71','71','72','72','72','73','73','73','74','74','74','75','75','75','76','76','76','77','77','77','78','78','78','79',
                '79','79','80','80','80','81','81','81','82','82','82','83','83','83','84','84','84','85','85','85','86','86','86','86'],
    'loss' : [9.21e+3,9.21e+3,9.21e+3,8.41e+3,8.41e+3,8.41e+3,8.33e+3,8.33e+3,8.33e+3,8.28e+3,8.28e+3,8.28e+3,8.24e+3,8.24e+3,8.24e+3,8.21e+3,8.21e+3,8.21e+3,8.18e+3,8.18e+3,8.18e+3,8.16e+3,8.16e+3,8.16e+3,
               8.15e+3,8.15e+3,8.15e+3,8.13e+3,8.13e+3,8.13e+3,8.12e+3,8.12e+3,8.12e+3,8.11e+3,8.11e+3,8.11e+3,8.11e+3,8.11e+3,8.11e+3,8.1e+3,8.1e+3,8.1e+3,8.09e+3,8.09e+3,8.09e+3,8.09e+3,8.09e+3,8.09e+3,
               8.09e+3,8.09e+3,8.09e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.08e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,8.07e+3,
               8.07e+3,8.07e+3,8.07e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,
               8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.06e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,
               8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3,8.05e+3 ,8.05e+3,8.05e+3,
               8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,
               8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,
               8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,
               8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.04e+3,8.03e+3,8.03e+3,8.03e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,
               8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3,8.02e+3]})

parse_training.drop_duplicates()


from matplotlib.ticker import MultipleLocator

fig, axes = plt.subplots(
    nrows=1, ncols=2,
    figsize=(6 , 2),
    constrained_layout=True
)
#fig, ax = plt.subplots(figsize=(5, 3))
axes[1].plot(parse_training.epoch, parse_training.loss)

# Set major ticks every 5 epochs
axes[1].xaxis.set_major_locator(MultipleLocator(25))

# Optional: Add minor ticks (not labeled) at each epoch
axes[1].xaxis.set_minor_locator(MultipleLocator(1))

axes[1].set_xlabel('Epoch')
axes[1].set_ylabel('Loss')
axes[1].set_title('Parse - Training Loss per epoch')

#fig, ax = plt.subplots(figsize=(5, 3))
axes[0].plot(l.epoch, l.loss)

# Set major ticks every 5 epochs
axes[0].xaxis.set_major_locator(MultipleLocator(25))

# Optional: Add minor ticks (not labeled) at each epoch
axes[0].xaxis.set_minor_locator(MultipleLocator(1))

axes[0].set_xlabel('Epoch')
axes[0].set_ylabel('Loss')
axes[0].set_title('10X - Training Loss per epoch')

fig.savefig('learning_curve.svg', dpi = 300)



#### Subtype marker analysis - Figure 4

genes = ['CD3D', 'CD3G', 'CD3E', 'CD247', 'CD4', 'CD8A', 'CD8B',
         'IL7R', 'IL2RG', 'IL2RA', 'CCR7', 'SELL',  'LEF1', 'TCF7',
         'CXCR3', 'TBX21', 'CCR4', 'GATA3', 'RORC', 'CXCR5', 'FOXP3',
         'GZMA', 'GZMB', 'GZMH', 'GZMK', 'GZMM', 'PRF1', 'GNLY',
        'TRAV1-2', 'SLC4A10', 'PDCD1', 'LAG3', 'B3GAT1', 'CX3CR1', 'ITGAE']

def gene_umaps(genes):
    fig = plt.figure(figsize=(2,4))
    gs = fig.add_gridspec(2, len(genes))

    max_gene = {}
    for gene in genes:
        max1 = tadata[:, gene].X.max()
        max2 = adata[:, gene].X.max()

        if max1 < max2:
            max_gene[gene] = max2
        else:
            max_gene[gene] = max1

    axes_top = []
    for i in range(len(genes)):
        ax = fig.add_subplot(gs[0,i])
        axes_top.append(ax)

    for i, gene in enumerate(genes):
        sc.pl.umap(tadata, color=gene, vmin=0, vmax=max_gene[gene], size=2, ax=axes_top[i], show=False, title=f'{gene} (10X)')

    axes_bottom = []
    for i in range(len(genes)):
        ax = fig.add_subplot(gs[1,i])
        axes_bottom.append(ax)

    for i, gene in enumerate(genes):
        sc.pl.umap(adata, color=gene, vmin=0, vmax=max_gene[gene], size=2, ax=axes_bottom[i], show=False, title=f'{gene} (Parse)')

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    return fig


for gene in genes:
    fig = gene_umaps([gene])
    fig.savefig(f'{gene}.png', dpi = 300)
    plt.close(fig)
    print(f'saved {gene}.png')


##### Supplementary Table 3

# Define the dataset names
dataset1_name = "10X"    # Formerly "tadata"
dataset2_name = "Parse"  # Formerly "adata"

# Function to get summary statistics for a dataset
def get_gene_stats(adata, genes_list, dataset_name):
    stats_list = []
    
    for gene in genes_list:
        # Extract expression values for this gene
        if gene in adata.var_names:
            expr = adata[:, gene].X.toarray().flatten() if scipy.sparse.issparse(adata.X) else adata[:, gene].X.flatten()
            
            # Calculate statistics
            stats = {
                'Gene': gene,
                'Dataset': dataset_name,
                'Mean': np.mean(expr),
                'Median': np.median(expr),
                'Std Dev': np.std(expr),
                'Min': np.min(expr),
                'Max': np.max(expr),
                '25th %ile': np.percentile(expr, 25),
                '75th %ile': np.percentile(expr, 75),
                'Nonzero %': (expr > 0).sum() / len(expr) * 100,
                'Cell Count': len(expr)
            }
            stats_list.append(stats)
        else:
            print(f"Warning: Gene {gene} not found in dataset")
    
    return pd.DataFrame(stats_list)

# Get statistics for both datasets
# Note: still using tadata and adata as the variable names, but with new labels
dataset1_stats = get_gene_stats(tadata, genes, dataset1_name)
dataset2_stats = get_gene_stats(adata, genes, dataset2_name)

# Combine into a single DataFrame
combined_stats = pd.concat([dataset1_stats, dataset2_stats])

# Sort by gene first, then by dataset
combined_stats = combined_stats.sort_values(['Gene', 'Dataset'])

# Reorder the columns for better readability
column_order = ['Gene', 'Dataset', 'Mean', 'Median', 'Std Dev', 'Min', 'Max', 
                '25th %ile', '75th %ile', 'Nonzero %', 'Cell Count']
combined_stats = combined_stats[column_order]

# Format the table
formatted_stats = combined_stats.round(3)
print(formatted_stats)

# Create a pivot table for easier visualization of gene-focused statistics
pivot_table = combined_stats.pivot(index='Gene', columns='Dataset', 
                                  values=['Mean', 'Median', 'Std Dev', 'Nonzero %'])

# Reorder columns for better readability
pivot_table = pivot_table.reindex(columns=[('Mean', dataset1_name), ('Mean', dataset2_name),
                                          ('Median', dataset1_name), ('Median', dataset2_name),
                                          ('Std Dev', dataset1_name), ('Std Dev', dataset2_name),
                                          ('Nonzero %', dataset1_name), ('Nonzero %', dataset2_name)])

print("\nPivot Table (Gene-focused):")
print(pivot_table.round(3))

# Save the pivot table to CSV
pivot_table.round(3).to_csv('gene_expression_pivot_table.csv')

# Create a gene-by-gene comparison table with the new dataset names
comparison_df = pd.DataFrame()

for gene in genes:
    dataset1_row = dataset1_stats[dataset1_stats['Gene'] == gene].iloc[0]
    dataset2_row = dataset2_stats[dataset2_stats['Gene'] == gene].iloc[0]
    
    # Create a row for this gene
    gene_comparison = {
        'Gene': gene,
        f'{dataset1_name} Mean': dataset1_row['Mean'],
        f'{dataset2_name} Mean': dataset2_row['Mean'],
        'Mean Ratio (10X/Parse)': dataset1_row['Mean'] / dataset2_row['Mean'] if dataset2_row['Mean'] != 0 else float('inf'),
        f'{dataset1_name} Median': dataset1_row['Median'],
        f'{dataset2_name} Median': dataset2_row['Median'],
        f'{dataset1_name} Std Dev': dataset1_row['Std Dev'],
        f'{dataset2_name} Std Dev': dataset2_row['Std Dev'],
        'Std Dev Ratio (10X/Parse)': dataset1_row['Std Dev'] / dataset2_row['Std Dev'] if dataset2_row['Std Dev'] != 0 else float('inf'),
        f'{dataset1_name} Nonzero %': dataset1_row['Nonzero %'],
        f'{dataset2_name} Nonzero %': dataset2_row['Nonzero %'],
        'Nonzero % Diff (10X-Parse)': dataset1_row['Nonzero %'] - dataset2_row['Nonzero %']
    }
    comparison_df = pd.concat([comparison_df, pd.DataFrame([gene_comparison])], ignore_index=True)

print("\nGene-by-Gene Comparison:")
print(comparison_df.round(3))

html_table = comparison_df.round(3).to_html()
with open('gene_expression_statistics.html', 'w') as f:
    f.write(html_table)

#### Supplementary Figure 2d

expression_diff = {
    'CD3D' : 5.654,
    'CD3G' : 1.791,
    'CD3E' : 3.017,
    'CD247' : 0.298,
    'CD4' : 0.594,
    'CD8A' : 2.400,
    'CD8B' : 3.107,
    'IL7R' : 0.809,
    'IL2RG' : 2.896,
    'IL2RA' : 0.253,
    'CCR7' : 0.713,
    'SELL' : 1.052,
    'LEF1' : 0.354,
    'TCF7' : 0.600,
    'CXCR3' : 14.711,
    'TBX21' : 1.560,
    'CCR4' : 0.701,
    'GATA3' : 1.016,
    'RORC' : 0.613,
    'CXCR5' : 0.761,
    'FOXP3' : 1.291,
    'GZMA' : 6.834,
    'GZMB' : 7.995,
    'GZMH' : 7.318,
    'GZMK' : 6.172,
    'GZMM' : 17.364,
    'PRF1' : 3.064,
    'GNLY' : 2.249,
    'TRAV1-2' : 17.078,
    'SLC4A10' : 0.255,
    'CX3CR1' : 1.471,
    'ITGAE' : 1.748,
    'PDCD1' : 1.693,
    'LAG3' : 8.414,
    'B3GAT1' : 0.123
}

def translate(name, ens):
    return ens.genes_by_name(name)[0].gene_id

# Create gene_id dictionary by translating gene names to Ensembl IDs
gene_id = {}
for gene_name in expression_diff.keys():
    try:
        gene_id[gene_name] = translate(gene_name, ens109)
    except:
        print(f"Could not translate {gene_name}")

with open('/Path/to/transcript_lengths_ensembl109.json', 'r') as f:
    lengths_109 = json.load(f)

# Get lengths for genes that exist in both dictionaries
lengths_for_genes = {}
for gene_name, ens_id in gene_id.items():
    if ens_id in lengths_109:
        lengths_for_genes[gene_name] = lengths_109[ens_id]

# Create data dictionary with matching genes
gene_names = list(lengths_for_genes.keys())
data = {
    'gene': gene_names,
    'expression_diff': [expression_diff[gene] for gene in gene_names],
    'length': [lengths_for_genes[gene] for gene in gene_names]
}

print(f"Found {len(data['length'])} genes with both expression and length data")
print(f"Expression diff range: {min(data['expression_diff']):.2f} - {max(data['expression_diff']):.2f}")
print(f"Length range: {min(data['length'])} - {max(data['length'])} bp")

# Pearson correlation (assumes linear relationship)
pearson_r, pearson_p = pearsonr(data['length'], data['expression_diff'])
print(f"\nPearson correlation: r = {pearson_r:.4f}, p-value = {pearson_p:.4f}")

# Spearman correlation (non-parametric, rank-based)
spearman_r, spearman_p = spearmanr(data['length'], data['expression_diff'])
print(f"Spearman correlation: rho = {spearman_r:.4f}, p-value = {spearman_p:.4f}")

# Create scatter plot with labels and colors
fig, ax = plt.subplots(figsize=(6, 5))

# Color code by expression difference
colors = ['navy' if expr_diff > 1 else 'firebrick' for expr_diff in data['expression_diff']]

ax.scatter(data['length'], data['expression_diff'], c=colors, alpha=1, s=50, edgecolors='black', linewidth=0.5)

# Add gene labels to each point
for i, gene in enumerate(data['gene']):
    ax.annotate(gene, 
                (data['length'][i], data['expression_diff'][i]),
                xytext=(3, 2),  # Offset the label slightly
                textcoords='offset points',
                fontsize=8,
                alpha=1)

ax.axhline(y=1, color='red', linestyle='--', linewidth=2, label='Expression diff = 1', alpha = 0.5)
ax.set_xlabel('Transcript Length (bp)')
ax.set_ylabel('Expression Difference')
ax.set_title(f'Spearman rho = {spearman_r:.4f}, p = {spearman_p:.4f}')
ax.grid(alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='navy', edgecolor='black', label='Higher in 10X (>1)'),
                   Patch(facecolor='firebrick', edgecolor='black', label='Higher in Parse (<1)')]
ax.legend(handles=legend_elements)

plt.tight_layout()
fig.savefig('subtype_genes_corr.svg', bbox_inches='tight', dpi=300)
plt.show()
