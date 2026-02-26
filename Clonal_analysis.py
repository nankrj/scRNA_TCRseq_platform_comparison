import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import ttest_rel, wilcoxon
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import euclidean


padata = sc.read_h5ad('/Path/to/PARSE_postQC_Tcells.h5ad')
tadata = sc.read_h5ad('/Path/to/10X_postQC_Tcells.h5ad')

import re
def extract_first_sequence(cdr3s_aa):
    # Regular expression to match the first sequence after 'TRB:' and before the first ';'
    match = re.search(r'TRB:([^;]+)', cdr3s_aa)
    if match:
        return match.group(1)  # Return the matched sequence without 'TRB:'
    else:
        return None  # In case no match is found


tadata.obs['TRB_cdr3_aa'] = tadata.obs.cdr3s_aa.apply(extract_first_sequence)

def extract_first_sequence(cdr3s_aa):
    # Regular expression to match the first sequence after 'TRB:' and before the first ';'
    match = re.search(r'TRA:([^;]+)', cdr3s_aa)
    if match:
        return match.group(1)  # Return the matched sequence without 'TRB:'
    else:
        return None  # In case no match is found

tadata.obs['TRA_cdr3_aa'] = tadata.obs.cdr3s_aa.apply(extract_first_sequence)

P1 = padata.obs[padata.obs.fig_label == '1'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P2 = padata.obs[padata.obs.fig_label == '2'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P3 = padata.obs[padata.obs.fig_label == '3'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P4 = padata.obs[padata.obs.fig_label == '4'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P5 = padata.obs[padata.obs.fig_label == '5'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P6 = padata.obs[padata.obs.fig_label == '6'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
P7 = padata.obs[padata.obs.fig_label == '7'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]

T1 = tadata.obs[tadata.obs.fig_label == '1'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T2 = tadata.obs[tadata.obs.fig_label == '2'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T3 = tadata.obs[tadata.obs.fig_label == '3'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T4 = tadata.obs[tadata.obs.fig_label == '4'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T5 = tadata.obs[tadata.obs.fig_label == '5'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T6 = tadata.obs[tadata.obs.fig_label == '6'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]
T7 = tadata.obs[tadata.obs.fig_label == '7'][['TRB_cdr3_aa', 'TRA_cdr3_aa']]

C_TRB_1 = set(P1.TRB_cdr3_aa.unique()) & set(T1.TRB_cdr3_aa.unique())
C_TRB_2 = set(P2.TRB_cdr3_aa.unique()) & set(T2.TRB_cdr3_aa.unique())
C_TRB_3 = set(P3.TRB_cdr3_aa.unique()) & set(T3.TRB_cdr3_aa.unique())
C_TRB_4 = set(P4.TRB_cdr3_aa.unique()) & set(T4.TRB_cdr3_aa.unique())
C_TRB_5 = set(P5.TRB_cdr3_aa.unique()) & set(T5.TRB_cdr3_aa.unique())
C_TRB_6 = set(P6.TRB_cdr3_aa.unique()) & set(T6.TRB_cdr3_aa.unique())
C_TRB_7 = set(P7.TRB_cdr3_aa.unique()) & set(T7.TRB_cdr3_aa.unique())

C_TRA_1 = set(P1.TRA_cdr3_aa.unique()) & set(T1.TRA_cdr3_aa.unique())
C_TRA_2 = set(P2.TRA_cdr3_aa.unique()) & set(T2.TRA_cdr3_aa.unique())
C_TRA_3 = set(P3.TRA_cdr3_aa.unique()) & set(T3.TRA_cdr3_aa.unique())
C_TRA_4 = set(P4.TRA_cdr3_aa.unique()) & set(T4.TRA_cdr3_aa.unique())
C_TRA_5 = set(P5.TRA_cdr3_aa.unique()) & set(T5.TRA_cdr3_aa.unique())
C_TRA_6 = set(P6.TRA_cdr3_aa.unique()) & set(T6.TRA_cdr3_aa.unique())
C_TRA_7 = set(P7.TRA_cdr3_aa.unique()) & set(T7.TRA_cdr3_aa.unique())

C_Pair_1 = set(zip(P1.TRB_cdr3_aa, P1.TRA_cdr3_aa)) & set(zip(T1.TRB_cdr3_aa, T1.TRA_cdr3_aa))
C_Pair_2 = set(zip(P2.TRB_cdr3_aa, P2.TRA_cdr3_aa)) & set(zip(T2.TRB_cdr3_aa, T2.TRA_cdr3_aa))
C_Pair_3 = set(zip(P3.TRB_cdr3_aa, P3.TRA_cdr3_aa)) & set(zip(T3.TRB_cdr3_aa, T3.TRA_cdr3_aa))
C_Pair_4 = set(zip(P4.TRB_cdr3_aa, P4.TRA_cdr3_aa)) & set(zip(T4.TRB_cdr3_aa, T4.TRA_cdr3_aa))
C_Pair_5 = set(zip(P5.TRB_cdr3_aa, P5.TRA_cdr3_aa)) & set(zip(T5.TRB_cdr3_aa, T5.TRA_cdr3_aa))
C_Pair_6 = set(zip(P6.TRB_cdr3_aa, P6.TRA_cdr3_aa)) & set(zip(T6.TRB_cdr3_aa, T6.TRA_cdr3_aa))
C_Pair_7 = set(zip(P7.TRB_cdr3_aa, P7.TRA_cdr3_aa)) & set(zip(T7.TRB_cdr3_aa, T7.TRA_cdr3_aa))

def common_clone(row, common_trb, common_tra, common_pairs):
    """
    For a given row (with TRB and TRA sequences), returns:
      - 'Pair' if the exact (TRB, TRA) pair exists in both datasets,
      - 'Both' if both sequences are common individually but not as a pair,
      - 'TRB' if only the TRB sequence is common,
      - 'TRA' if only the TRA sequence is common,
      - 'None' if neither is common.
    """
    trb = row['TRB_cdr3_aa']
    tra = row['TRA_cdr3_aa']
    
    # Check if the exact pair is common.
    if (trb, tra) in common_pairs:
        return 'Pair'
    # Check if both chains are common individually.
    elif (trb in common_trb) and (tra in common_tra):
        return 'Partly'
    elif trb in common_trb:
        return 'Partly'
    elif tra in common_tra:
        return 'Partly'
    else:
        return 'None'


P1['common'] = P1.apply(lambda row: common_clone(row, C_TRB_1, C_TRA_1, C_Pair_1), axis=1)
T1['common'] = T1.apply(lambda row: common_clone(row, C_TRB_1, C_TRA_1, C_Pair_1), axis=1)
P2['common'] = P2.apply(lambda row: common_clone(row, C_TRB_2, C_TRA_2, C_Pair_2), axis=1)
T2['common'] = T2.apply(lambda row: common_clone(row, C_TRB_2, C_TRA_2, C_Pair_2), axis=1)
P3['common'] = P3.apply(lambda row: common_clone(row, C_TRB_3, C_TRA_3, C_Pair_3), axis=1)
T3['common'] = T3.apply(lambda row: common_clone(row, C_TRB_3, C_TRA_3, C_Pair_3), axis=1)
P4['common'] = P4.apply(lambda row: common_clone(row, C_TRB_4, C_TRA_4, C_Pair_4), axis=1)
T4['common'] = T4.apply(lambda row: common_clone(row, C_TRB_4, C_TRA_4, C_Pair_4), axis=1)
P5['common'] = P5.apply(lambda row: common_clone(row, C_TRB_5, C_TRA_5, C_Pair_5), axis=1)
T5['common'] = T5.apply(lambda row: common_clone(row, C_TRB_5, C_TRA_5, C_Pair_5), axis=1)
P6['common'] = P6.apply(lambda row: common_clone(row, C_TRB_6, C_TRA_6, C_Pair_6), axis=1)
T6['common'] = T6.apply(lambda row: common_clone(row, C_TRB_6, C_TRA_6, C_Pair_6), axis=1)
P7['common'] = P7.apply(lambda row: common_clone(row, C_TRB_7, C_TRA_7, C_Pair_7), axis=1)
T7['common'] = T7.apply(lambda row: common_clone(row, C_TRB_7, C_TRA_7, C_Pair_7), axis=1)

p1_cell_count = P1.value_counts().reset_index()
p2_cell_count = P2.value_counts().reset_index()
p3_cell_count = P3.value_counts().reset_index()
p4_cell_count = P4.value_counts().reset_index()
p5_cell_count = P5.value_counts().reset_index()
p6_cell_count = P6.value_counts().reset_index()
p7_cell_count = P7.value_counts().reset_index()
t1_cell_count = T1.value_counts().reset_index()
t2_cell_count = T2.value_counts().reset_index()
t3_cell_count = T3.value_counts().reset_index()
t4_cell_count = T4.value_counts().reset_index()
t5_cell_count = T5.value_counts().reset_index()
t6_cell_count = T6.value_counts().reset_index()
t7_cell_count = T7.value_counts().reset_index()

#Saving the counts to make bubble plots in R

p1_cell_count.to_csv('p1_clone_counts.csv')
p2_cell_count.to_csv('p2_clone_counts.csv')
p3_cell_count.to_csv('p3_clone_counts.csv')
p4_cell_count.to_csv('p4_clone_counts.csv')
p5_cell_count.to_csv('p5_clone_counts.csv')
p6_cell_count.to_csv('p6_clone_counts.csv')
p7_cell_count.to_csv('p7_clone_counts.csv')
t1_cell_count.to_csv('t1_clone_counts.csv')
t2_cell_count.to_csv('t2_clone_counts.csv')
t3_cell_count.to_csv('t3_clone_counts.csv')
t4_cell_count.to_csv('t4_clone_counts.csv')
t5_cell_count.to_csv('t5_clone_counts.csv')
t6_cell_count.to_csv('t6_clone_counts.csv')
t7_cell_count.to_csv('t7_clone_counts.csv')


com_count_P1 = P1.common.value_counts().rename('P1')
com_count_P2 = P2.common.value_counts().rename('P2')
com_count_P3 = P3.common.value_counts().rename('P3')
com_count_P4 = P4.common.value_counts().rename('P4')
com_count_P5 = P5.common.value_counts().rename('P5')
com_count_P6 = P6.common.value_counts().rename('P6')
com_count_P7 = P7.common.value_counts().rename('P7')
com_count_T1 = T1.common.value_counts().rename('T1')
com_count_T2 = T2.common.value_counts().rename('T2')
com_count_T3 = T3.common.value_counts().rename('T3')
com_count_T4 = T4.common.value_counts().rename('T4')
com_count_T5 = T5.common.value_counts().rename('T5')
com_count_T6 = T6.common.value_counts().rename('T6')
com_count_T7 = T7.common.value_counts().rename('T7')

dfs = [com_count_T1, com_count_T2, com_count_T3, com_count_T4, com_count_T5, com_count_T6, com_count_T7, com_count_P1, com_count_P2, com_count_P3, com_count_P4, com_count_P5, com_count_P6, com_count_P7]


# Merge all counts
combined_df = pd.concat(dfs, axis=1).fillna(0).astype(int)

combined_df.T

sns.set_theme(style="white")

categories = [
    'Pair',
    'Partly',
    'None'
]

# Convert to percentage
combined_percent = combined_df.div(combined_df.sum(axis=0), axis=1) * 100

# Reorder rows based on biotype order
combined_percent = combined_percent.loc[categories]
combined_df = combined_df.loc[categories]

# Define colors for each biotype
biotype_colors = {
    'Pair' : '#0A3431',
    'Partly' : '#2A6858',
    'None' : '#4A9B7F'
}

#### Figure 3e,f)

# Function to create stacked barplots with seaborn styling
def plot_stacked_barplot(df, color_dict, legend_labels=None, figsize=(12, 6), width=0.85, edgecolor='black', alpha=0.8):

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    
    # Get all biotypes and categories
    biotypes = df.index
    categories = df.columns
    
    # Get color list for biotypes
    colors = [color_dict.get(biotype, "#999999") for biotype in biotypes]
    
    # Initialize bottom positions for stacking
    bottoms = np.zeros(len(categories))
    
    # Plot each biotype as a segment of the stacked bar
    for i, biotype in enumerate(biotypes):
        values = df.loc[biotype].values
        # Plot the bars
        ax.bar(categories, values, bottom=bottoms, width=width, 
               label=biotype, color=colors[i], edgecolor=edgecolor, alpha=alpha)
        # Update the bottoms for the next biotype
        bottoms += values
    
    # Identify spacer columns
    spacer_cols = [col for col in categories if "Spacer" in col]
    
    # Define main categories and their positions
    main_categories = ["1", '2', '3', '4', '5', '6', '7']
    
    # Find positions for each main category group
    group_centers = []
    current_group = 0
    group_start_idx = 0
    
    for i, col in enumerate(categories):
        # Check if we've reached a spacer or the end
        is_last_column = (i == len(categories) - 1)
        is_spacer = "Spacer" in col
        
        if is_spacer or is_last_column:
            # Calculate group center before moving to next group
            if is_last_column and not is_spacer:
                end_idx = i  # Include this column in the group
            else:
                end_idx = i - 1  # Don't include the spacer
                
            if end_idx >= group_start_idx:  # Make sure we have at least one column
                center = (group_start_idx + end_idx) / 2
                group_centers.append(center)
                
            if is_spacer and current_group < len(main_categories) - 1:
                current_group += 1
                group_start_idx = i + 1  # Start next group after spacer
    
    # If we have fewer centers than categories, add the last one
    if len(group_centers) < len(main_categories):
        remaining_categories = len(main_categories) - len(group_centers)
        for i in range(remaining_categories):
            # Calculate position for the last group if missed
            last_spacer_idx = max([i for i, col in enumerate(categories) if "Spacer" in col], default=-1)
            if last_spacer_idx < len(categories) - 1:
                center = (last_spacer_idx + 1 + len(categories) - 1) / 2
                group_centers.append(center)

    # Format the x-tick labels (remove "all_", ">20 cells_", etc. prefixes)
    clean_labels = []
    for col in categories:
        if "Spacer" in col:
            clean_labels.append(" ")
        else:
            # Extract the subcategory part (Common, PARSE, 10X)
            if "T" in col:
                clean_labels.append("10X")
            elif "P" in col:
                clean_labels.append("PARSE")
            else:
                clean_labels.append(col)  # Default fallback
    
    # Set x-tick positions and labels
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels(clean_labels, rotation = 45, ha='right')
    
    # Add main category labels above
    for i, (pos, label) in enumerate(zip(group_centers, main_categories)):
        ax.text(pos, -0.13, label, ha='center', va='top', fontsize=12, 
                fontweight='bold', transform=ax.get_xaxis_transform())
    
    # Clean up legend labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Apply custom legend labels if provided
    if legend_labels:
        clean_legend_labels = [legend_labels.get(label, label) for label in labels]
    else:
        # Default: replace underscores with spaces and capitalize first letter
        clean_legend_labels = [label.replace('_', ' ').capitalize() for label in labels]
    
    # Formatting
    ax.set_ylabel("Absolute Counts")
    #ax.set_title("Biotype Distribution")
    ax.legend(handles, clean_legend_labels, title="", bbox_to_anchor=(1.05, 1), loc="upper left")
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    
    return fig, ax


for space_col in ["Spacer_1", "Spacer_2", "Spacer_3", "Spacer_4", "Spacer_5", "Spacer_6"]:
    combined_percent[space_col] = float("nan")  # Fill with NaN
    combined_df[space_col] = float('nan')

# Define the order of columns with explicit NaN spacer columns
ordered_columns = [col for col in combined_percent.columns if "1" in col] + \
                  [col for col in combined_percent.columns if "2" in col] + \
                  [col for col in combined_percent.columns if "3" in col] + \
                  [col for col in combined_percent.columns if "4" in col] + \
                  [col for col in combined_percent.columns if "5" in col] + \
                  [col for col in combined_percent.columns if "6" in col] + \
                  [col for col in combined_percent.columns if "7" in col]

# Reorder the dataframe based on the new column order
combined_percent = combined_percent[ordered_columns]
combined_df = combined_df[ordered_columns]

# Handle missing colors by assigning a default gray color
color_list = [biotype_colors.get(b, "#999999") for b in combined_percent.index]

legend_labels = {
    'Pair' : 'Overlap',
    'Partly' : 'Partly overlap',
    'None' : 'Unique'}

fig, ax = plot_stacked_barplot(combined_df, biotype_colors, legend_labels = legend_labels)

plt.savefig('Common_barplot.svg', bbox_inches='tight', dpi=300)

fig, ax = plot_stacked_barplot(combined_percent, biotype_colors, legend_labels = legend_labels)

plt.savefig('Common_barplot_percent.svg', bbox_inches='tight', dpi=300)


##### Calculating diversity
def calc_shannon_norm(df):
    df_copy = df.copy()
    df_copy['frac'] = df_copy['count'] / df_copy['count'].sum()
    # Handle log(0) by replacing zeros with a small value, or excluding them
    nonzero_frac = df_copy['frac'][df_copy['frac'] > 0]
    return -sum(nonzero_frac * np.log(nonzero_frac))


def calc_d50(df):
    df_copy = df.copy().sort_values('count', ascending=False)
    df_copy['cumfrac'] = df_copy['count'].cumsum() / df_copy['count'].sum()
    
    # Count clones needed for 50%
    n50 = (df_copy['cumfrac'] < 0.5).sum() + 1
    S = len(df_copy)
    
    return (n50 / S) * 100

def calc_pielou(df, base=None):
    if base is None:
        base = len(df)  # Default normalization by number of species
    shannon = calc_shannon_norm(df)
    return shannon / np.log(base)

def calc_simpson_diversity(df):
    N = df['count'].sum()
    numerator = np.sum(df['count'] * (df['count'] - 1))
    denominator = N * (N - 1)
    return 1 - (numerator / denominator) if denominator > 0 else 0

def calc_chao1(df):
    S = len(df)
    singletons = sum(df['count'] == 1)
    doubletons = sum(df['count'] == 2)
    if doubletons == 0:
        return S + singletons  # Avoid division by zero
    return S + (singletons ** 2) / (2 * doubletons)

def calc_gini(df):
    values = np.sort(df['count'].values)
    n = len(values)
    cum_vals = np.cumsum(values)
    return (2 * np.sum((np.arange(1, n+1) * values))) / (n * np.sum(values)) - (n + 1) / n

# Define sample pairs (10X and PARSE samples)
sample_pairs = [
    ('T1', 'P1', t1_cell_count, p1_cell_count),
    ('T2', 'P2', t2_cell_count, p2_cell_count),
    ('T3', 'P3', t3_cell_count, p3_cell_count),
    ('T4', 'P4', t4_cell_count, p4_cell_count),
    ('T5', 'P5', t5_cell_count, p5_cell_count),
    ('T6', 'P6', t6_cell_count, p6_cell_count),
    ('T7', 'P7', t7_cell_count, p7_cell_count),
]

# Define diversity metrics functions
diversity_metrics = {
    'shannon': calc_shannon_norm,
    'pielou': calc_pielou,
    'simpson_d': calc_simpson_diversity,
    'chao1': calc_chao1,
    'gini': calc_gini,
    'd50': calc_d50
}

# Initialize results dictionary with sample names
results = {'sample': []}
for metric_name in diversity_metrics:
    results[metric_name] = []

# Calculate all metrics for each sample
for t_name, p_name, t_data, p_data in sample_pairs:
    # Add sample names
    results['sample'].extend([t_name, p_name])
    
    # Calculate each metric and store
    for metric_name, metric_func in diversity_metrics.items():
        # Make copies to avoid modifying original dataframes
        t_copy = t_data.copy() 
        p_copy = p_data.copy()
        
        # Calculate and store results
        results[metric_name].append(metric_func(t_copy))
        results[metric_name].append(metric_func(p_copy))

# Create the final dataframe
div_df = pd.DataFrame(results)
print(div_df)

#### Figure 3g-l)

# Create a new column to separate Tumor and Paired samples
div_df['Type'] = ['10X' if 'T' in s else 'PARSE' for s in div_df['sample']]

# Convert 'sample' to numeric order for proper plotting
div_df['Sample Number'] = div_df['sample'].str.extract(r'(\d+)').astype(int)

# Define metrics to visualize
metrics = [
    ('shannon', 'Shannon', 'Shannon Diversity Index', 5, 8.6),
    ('pielou', 'Pielou', 'Pielou Evenness Index', 0.8, 1.01),
    ('simpson_d', 'Simpson', 'Simpsons Diversity Index', 0.95, 1.0035),
    ('gini', 'Gini', 'Gini Index', None, 1.05),
    ('chao1', 'Chao1', 'Chao1 Species Richness Estimator', 0, 70000),
    ('d50', 'D50', 'D50', None, None)
]

# Plot
fig, axes = plt.subplots(
    nrows=2, ncols=6,
    figsize=(16, 4),
    constrained_layout=True,
    gridspec_kw={'height_ratios': [1, 1.5]}
)

# Flatten for easier iteration
#axes = axes.flatten()

# Create paired plots for each metric (zoomed and full range)
for i, (metric, short_label, full_label, min_val, max_val) in enumerate(metrics):
    # Zoomed view (left column)
    left_idx = [0,i]
    sns.barplot(x='Sample Number', y=metric, hue='Type', data=div_df, 
                palette={'10X': 'navy', 'PARSE': 'firebrick'}, 
                edgecolor='black', alpha=0.8, ax=axes[0,i])
    axes[0,i].set_xlabel('')
    axes[0,i].set_ylabel(short_label)
    axes[0,i].grid(False)
    if min_val is not None and max_val is not None:
        axes[0,i].set_ylim(min_val, max_val)
    if i == 5:
        axes[0,i].legend_.remove()#legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        axes[0,i].legend_.remove()
    
    # Full range view (right column)
    right_idx = [1,i]
    sns.barplot(x='Sample Number', y=metric, hue='Type', data=div_df, 
                palette={'10X': 'navy', 'PARSE': 'firebrick'}, 
                edgecolor='black', alpha=0.8, ax=axes[1,i])
    axes[1,i].set_xlabel('')
    axes[1,i].set_ylabel(full_label)
    axes[1,i].grid(False)
    if min_val is None and max_val is not None:
        axes[1,i].set_ylim(min_val, max_val)
    axes[1,i].legend_.remove()


plt.tight_layout()
plt.savefig('Diversity_barplot.svg', bbox_inches='tight', dpi=300)

#### Correlation and agreement
metrics = ['shannon', 'pielou', 'simpson_d', 'gini', 'chao1', 'd50']

def create_paired_data(df, metrics):
    """Create paired dataframe for method comparison"""
    paired_data = []
    
    for patient in df['Sample Number'].unique():
        patient_data = df[df['Sample Number'] == patient]
        if len(patient_data) == 2 and set(patient_data['Type']) == {'10X', 'PARSE'}:
            tenx_data = patient_data[patient_data['Type'] == '10X'][metrics].iloc[0]
            parse_data = patient_data[patient_data['Type'] == 'PARSE'][metrics].iloc[0]
            
            for metric in metrics:
                paired_data.append({
                    'Patient': patient,
                    'Metric': metric,
                    '10X': tenx_data[metric],
                    'PARSE': parse_data[metric],
                    'Difference': parse_data[metric] - tenx_data[metric],
                    'Average': (tenx_data[metric] + parse_data[metric]) / 2,
                    'Relative_Difference': (parse_data[metric] - tenx_data[metric]) / tenx_data[metric] * 100 if tenx_data[metric] != 0 else np.nan
                })
    
    return pd.DataFrame(paired_data)

paired_df = create_paired_data(div_df, metrics)

correlation_results = []
for metric in metrics:
    metric_data = paired_df[paired_df['Metric'] == metric]
    if len(metric_data) > 2:
        # Calculate correlations
        pearson_r, pearson_p = stats.pearsonr(metric_data['10X'], metric_data['PARSE'])
        spearman_r, spearman_p = stats.spearmanr(metric_data['10X'], metric_data['PARSE'])
        
        correlation_results.append({
            'Metric': metric,
            'Pearson_r': pearson_r,
            'Pearson_p': pearson_p,
            'Spearman_r': spearman_r,
            'Spearman_p': spearman_p
        })
        
        print(f"{metric:12}: Pearson r={pearson_r:.3f} (p={pearson_p:.3f}), "
              f"Spearman r={spearman_r:.3f} (p={spearman_p:.3f})")

def calculate_icc(data):
    """Calculate ICC(2,1) - two-way random effects, single measures, absolute agreement"""
    n_subjects = len(data)
    n_raters = 2
    
    # Reshape data for ICC calculation
    ratings = np.array([data['10X'].values, data['PARSE'].values]).T
    
    # Calculate mean squares
    total_mean = np.mean(ratings)
    
    # Between subjects sum of squares
    subject_means = np.mean(ratings, axis=1)
    MSB = n_raters * np.sum((subject_means - total_mean) ** 2) / (n_subjects - 1)
    
    # Within subjects sum of squares  
    MSW = np.sum((ratings - subject_means.reshape(-1, 1)) ** 2) / (n_subjects * (n_raters - 1))
    
    # Between raters sum of squares
    rater_means = np.mean(ratings, axis=0)
    MSR = n_subjects * np.sum((rater_means - total_mean) ** 2) / (n_raters - 1)
    
    # Error sum of squares
    MSE = (np.sum((ratings - subject_means.reshape(-1, 1) - rater_means + total_mean) ** 2) / 
           ((n_subjects - 1) * (n_raters - 1)))
    
    # ICC calculation
    icc = (MSB - MSE) / (MSB + (n_raters - 1) * MSE)
    
    return max(0, min(1, icc))  # Constrain between 0 and 1

# Calculate ICC for each metric and add to correlation_results
for i, metric in enumerate(metrics):
    metric_data = paired_df[paired_df['Metric'] == metric]
    if len(metric_data) > 2:
        icc_value = calculate_icc(metric_data)
        correlation_results[i]['ICC'] = icc_value
        print(f"{metric:12}: ICC = {icc_value:.3f}")

# Create scaled version for distance calculations only
div_df_scaled = div_df.copy()
scaler = MinMaxScaler()
div_df_scaled[metrics] = scaler.fit_transform(div_df[metrics])

print("Data scaled using MinMaxScaler for distance calculations...")

# Calculate euclidean distances and cosine similarities on scaled data
euclidean_distances = []
cosine_similarities = []
patients_with_both_methods = []

for patient, group in div_df_scaled.groupby('Sample Number'):
    if group.shape[0] == 2:
        group = group.sort_values(by='Type')
        if '10X' in group['Type'].values and 'PARSE' in group['Type'].values:
            vec_10X = group[group['Type'] == '10X'][metrics].values[0]
            vec_Parse = group[group['Type'] == 'PARSE'][metrics].values[0]

            euc_dist = euclidean(vec_10X, vec_Parse)
            cos_sim = cosine_similarity(vec_10X.reshape(1, -1), vec_Parse.reshape(1, -1))[0][0]
            
            euclidean_distances.append(euc_dist)
            cosine_similarities.append(cos_sim)
            patients_with_both_methods.append(patient)
    else:
        print(f"Skipping patient {patient}: incomplete data.")

# Create results DataFrame
results_df = pd.DataFrame({
    'Patient': patients_with_both_methods,
    'Euclidean_Distance': euclidean_distances,
    'Euclidean_Agreement': [1/(1+euc) for euc in euclidean_distances],
    'Cosine_Similarity': cosine_similarities,
    'Cosine_Distance': [1 - sim for sim in cosine_similarities]
})

plt.style.use('default')
sns.set_style("whitegrid")
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 8
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['legend.fontsize'] = 7

#### Figure 3m,n)

def create_agreement_heatmap():
    """Create method agreement heatmap based on raw data correlations"""
    correlation_df = pd.DataFrame(correlation_results)
    
    fig, ax = plt.subplots(1, 1, figsize=(6, 3))
    
    # Prepare data for heatmap
    heatmap_data = correlation_df[['Metric', 'Pearson_r', 'Spearman_r', 'ICC']].set_index('Metric')
    
    # Create heatmap
    sns.heatmap(heatmap_data.T, annot=True, fmt='.3f', cmap='RdYlBu_r', 
                center=0.5, vmin=0, vmax=1, ax=ax, 
                cbar_kws={'label': 'Agreement Score'})
    
    ax.set_title('Method Agreement (Raw Values)', fontweight='bold')
    ax.set_yticklabels(['Pearson r', 'Spearman r', 'ICC'], rotation=0)
    
    plt.tight_layout()
    return fig

def create_distance_heatmap():
    """Create patient-level distance heatmap based on scaled data"""
    fig, ax = plt.subplots(1, 1, figsize=(6, 2.5))
    
    heatmap_data = results_df[['Patient', 'Euclidean_Agreement', 'Cosine_Similarity']].set_index('Patient')
    
    # Create heatmap
    sns.heatmap(heatmap_data.T, annot=True, fmt='.3f', cmap='RdYlBu_r', 
                center=0.5, vmin=0, vmax=1, ax=ax, 
                cbar_kws={'label': 'Agreement Score'})
    
    ax.set_title('Patient-Level Distance Agreement (Scaled Values)', fontweight='bold')
    ax.set_yticklabels(['Euclidean Agreement', 'Cosine Similarity'], rotation=0)
    
    plt.tight_layout()
    return fig

fig1 = create_agreement_heatmap()
fig2 = create_distance_heatmap()

fig1.savefig('method_agreement_heatmap.svg', dpi=300, bbox_inches='tight')
fig2.savefig('patient_distance_heatmap.svg', dpi=300, bbox_inches='tight')


#### Supplementary figure 2c
chao1_results = {f'P{i}': [] for i in range(1, 8)}
shannon_results = {f'P{i}': [] for i in range(1, 8)}

for i in range(10):
    adata_downsampled = padata.copy()
    sc.pp.subsample(adata_downsampled, n_obs=25479, random_state=i)
    
    for patient_num in range(1, 8):
        patient_data = adata_downsampled.obs[adata_downsampled.obs.fig_label == str(patient_num)][['TRB_cdr3_aa', 'TRA_cdr3_aa']].value_counts().reset_index()
        
        chao1 = calc_chao1(patient_data)
        shannon = calc_shannon_norm(patient_data)
        
        chao1_results[f'P{patient_num}'].append(chao1)
        shannon_results[f'P{patient_num}'].append(shannon)

# Calculate means and SDs
print("Chao1 results:")
for patient in sorted(chao1_results.keys()):
    mean_val = np.mean(chao1_results[patient])
    sd_val = np.std(chao1_results[patient])
    print(f"{patient}: {mean_val:.2f} ± {sd_val:.2f}")

print("\nShannon results:")
for patient in sorted(shannon_results.keys()):
    mean_val = np.mean(shannon_results[patient])
    sd_val = np.std(shannon_results[patient])
    print(f"{patient}: {mean_val:.2f} ± {sd_val:.2f}")

fig, axes = plt.subplots(
    nrows=1, ncols=2,
    figsize=(12, 5),
    constrained_layout=True
)

patients = ['1', '2', '3', '4', '5', '6', '7']

# Chao1 data
x10x_chao1 = [17157.50, 35594.17, 12519.34, 1543.47, 23077.57, 2927.15, 26938.60]
parse_orig_chao1 = [16699.40, 18764.36, 20516.11, 66013.48, 48023.01, 27386.82, 218923.00]
parse_down_chao1 = [15579.88, 15241.43, 14233.07, 61045.62, 44191.08, 21592.77, 189405.47]
parse_down_chao1_sd = [1638.73, 3694.88, 1515.62, 6124.07, 6276.16, 1448.41, 19137.15]

# Shannon data
x10x_shannon = [6.74, 6.30, 6.16, 6.22, 6.51, 5.18, 6.25]
parse_orig_shannon = [6.60, 6.17, 6.41, 7.81, 7.08, 6.58, 8.42]
parse_down_shannon = [6.41, 6.00, 6.24, 7.58, 6.89, 6.42, 8.18]
parse_down_shannon_sd = [0.03, 0.02, 0.02, 0.01, 0.02, 0.02, 0.01]

x = np.arange(len(patients))
width = 0.25

# Chao1
axes[0].bar(x - width, x10x_chao1, width, label='10X', alpha=0.8, color='navy', edgecolor='black')
axes[0].bar(x, parse_orig_chao1, width, label='Parse (Original)', alpha=0.8, color='firebrick', edgecolor='black')
axes[0].bar(x + width, parse_down_chao1, width, label='Parse (Downsampled)', 
           yerr=parse_down_chao1_sd, capsize=3, alpha=0.8, color='orange', edgecolor='black')
axes[0].set_title('Chao1', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Patient', fontsize=11)
axes[0].set_ylabel('Chao1 Index', fontsize=11)
axes[0].set_xticks(x)
axes[0].set_xticklabels(patients)
axes[0].grid(axis='y', alpha=0.3)
axes[0].legend(fontsize=10)

# Shannon
axes[1].bar(x - width, x10x_shannon, width, label='10X', alpha=0.8, color='navy', edgecolor='black')
axes[1].bar(x, parse_orig_shannon, width, label='Parse (Original)', alpha=0.8, color='firebrick', edgecolor='black')
axes[1].bar(x + width, parse_down_shannon, width, label='Parse (Downsampled)', 
           yerr=parse_down_shannon_sd, capsize=3, alpha=0.8, color='orange', edgecolor='black')
axes[1].set_title('Shannon', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Patient', fontsize=11)
axes[1].set_ylabel('Shannon Index', fontsize=11)
axes[1].set_xticks(x)
axes[1].set_xticklabels(patients)
axes[1].set_ylim(5, 8.5)
axes[1].grid(axis='y', alpha=0.3)
axes[1].legend(fontsize=10)

plt.tight_layout()
plt.savefig('supp_figure_downsampling_validation.svg', dpi=300, bbox_inches='tight')
plt.show()

