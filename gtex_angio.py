from scipy.stats import spearmanr, hypergeom
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from statsmodels.stats.multitest import multipletests
import textwrap

# ================== Angiogenesis marker list ==================
ANGIOGENESIS_MARKERS = [
    "VEGFA", "KDR", "ANGPT1", "PDGFB", "TIE1", "FLT1",
    "PECAM1", "CDH5", "NOS3", "ENG", "HIF1A", "CXCL8",
    "MMP2", "THBS1", "TGFB1", "EFNA1", "NRP1", "SPP1", "VEGFB",
    "FGF2", "FGF1", "PDGFA", "ANGPT2", "ANGPTL4", "DLL4",
    "NOTCH1", "NOTCH4", "EPHB4", "EFNB2", "COL18A1", "TIMP1",
    "TIMP3", "ICAM1", "VCAM1", "ADAMTS1", "ADAMTS9", "CXCR4",
    "CXCL12", "MMP9", "MMP14", "POSTN", "LOX", "PLAUR", "VEGFC",
    "VEGFD", "NRP2", "TIE2", "JAG1", "BMP4", "IL6"
]
# ===============================================================

def load_msigdb_gmt(gmt_file):
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            geneset_name = parts[0]
            genes = set(parts[2:])
            gene_sets[geneset_name] = genes
    return gene_sets

def process_gtex_file(file_path, target_gene_id, gene_id_to_name, control_marker_gene=None):
    """
    Processes a GTEx file and returns:
      - A DataFrame with correlation results for angiogenesis markers
      - A DataFrame with CPD and selected marker expression data (PECAM1 and control marker if available)
    """
    print(f"\nProcessing file: {os.path.basename(file_path)}")
    
    try:
        # Load data
        data = pd.read_csv(file_path, sep="\t", skiprows=2, index_col=0)
        # Create or update gene_id_to_name mapping if "Description" is available
        if "Description" in data.columns:
            gene_id_to_name = data["Description"].to_dict()
        data = data.drop(columns=["Description"], errors='ignore').apply(pd.to_numeric, errors="coerce").T
        data = data.loc[:, data.std() > 0]  # Remove constant genes

        # Check if CPD exists
        if target_gene_id not in data.columns:
            print(f"CPD ({target_gene_id}) not found")
            return None, None

        # Get angiogenesis markers present in dataset using gene names
        marker_ids = [g for g in data.columns if gene_id_to_name.get(g, g) in ANGIOGENESIS_MARKERS]
        if not marker_ids:
            print("No angiogenesis markers found in this dataset")
            return None, None
            
        # Calculate correlations (Spearman) for angiogenesis markers
        cpd_expression = data[target_gene_id]
        correlations = []
        for gene_id in marker_ids:
            marker_exp = data[gene_id]
            rho, pval = spearmanr(cpd_expression, marker_exp)
            correlations.append({
                'gene_id': gene_id,
                'gene_name': gene_id_to_name.get(gene_id, gene_id),
                'correlation': rho,
                'p_value': pval
            })

        # Create DataFrame and adjust p-values
        df = pd.DataFrame(correlations)
        df['fdr'] = multipletests(df['p_value'], method='fdr_bh')[1]

        # Build a DataFrame for scatter plots: include CPD and markers of interest
        marker_columns = {}
        # Try to extract PECAM1 expression if available
        for gene in data.columns:
            if gene_id_to_name.get(gene, gene) == 'PECAM1':
                marker_columns['PECAM1'] = gene
                break
        # Try to extract the control marker gene if specified
        if control_marker_gene:
            for gene in data.columns:
                if gene_id_to_name.get(gene, gene) == control_marker_gene:
                    marker_columns[control_marker_gene] = gene
                    break
        
        markers_df = None
        if marker_columns:
            markers_df = data[[target_gene_id] + list(marker_columns.values())].copy()
            # Rename columns: CPD for target_gene_id and use gene names for markers
            rename_dict = {target_gene_id: 'CPD'}
            for name, gene_id in marker_columns.items():
                rename_dict[gene_id] = name
            markers_df.rename(columns=rename_dict, inplace=True)
        
        return df.sort_values('correlation', ascending=False), markers_df

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None, None

def plot_angiogenesis_correlations(results_df, tissue_name):
    """Visualize correlations with angiogenesis markers"""
    plt.figure(figsize=(10, 6))
    
    # Sort values for plotting
    df = results_df.sort_values('correlation')
    
    # Color-code by significance
    colors = ['#1f77b4' if fdr < 0.05 else 'gray' for fdr in df['fdr']]
    
    plt.barh(
        y=df['gene_name'],
        width=df['correlation'],
        color=colors,
        edgecolor='black'
    )
    
    # Mark significant bars with '*'
    for i, (_, row) in enumerate(df.iterrows()):
        if row['fdr'] < 0.05:
            plt.text(row['correlation']/2, i, '*', ha='center', va='center', color='white', fontsize=12)
    
    plt.title(f"CPD Correlation with Angiogenesis Markers\n({tissue_name})")
    plt.xlabel("Spearman Correlation Coefficient")
    plt.ylabel("Gene")
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

def plot_scatter_markers(markers_df, tissue_name):
    """
    Visualize CPD vs marker expression for each selected marker (e.g., PECAM1 and control marker)
    in a single scatter plot with regression lines and correlation coefficients.
    """
    if markers_df is None or markers_df.empty:
        print("No marker data available for plotting.")
        return
    
    plt.figure(figsize=(8, 6))
    # Log-transform TPM values (add 1 to avoid log(0))
    data_log = np.log10(markers_df + 1)
    
    # Define colors for different markers
    marker_colors = ['blue', 'green', 'purple', 'orange']
    color_idx = 0
    for col in data_log.columns:
        if col == 'CPD':
            continue
        rho, pval = spearmanr(data_log['CPD'], data_log[col])
        sns.regplot(
            x='CPD', 
            y=col, 
            data=data_log,
            scatter_kws={'alpha': 0.5, 'color': marker_colors[color_idx % len(marker_colors)]},
            line_kws={'color': marker_colors[color_idx % len(marker_colors)], 'linestyle': '--'},
            label=f'{col} (r={rho:.2f}, p={pval:.2e})'
        )
        color_idx += 1

    plt.title(f'CPD vs Marker Expression in {tissue_name}\n(log10(TPM + 1))')
    plt.xlabel('CPD Expression (log10(TPM + 1))')
    plt.ylabel('Marker Expression (log10(TPM + 1))')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()

# Configuration
FILES = [
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_adrenal_gland.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_brain_spinal_cord_cervical_c-1.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_esophagus_muscularis.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_lung.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_whole_blood.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_cells_cultured_fibroblasts.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_minor_salivary_gland.gct.gz",
    "/Users/atharva/Downloads/tissue_tpm_data/gene_tpm_v10_thyroid.gct.gz",


    # ... other files ...
]

CPD_ID = "ENSG00000108582.12"  # Verify this is the correct CPD ID
CONTROL_MARKER = "MC2R"  # Change this to your marker gene of choice (e.g., VIM for vimentin)

if __name__ == "__main__":
    all_results = {}
    
    for file_path in FILES:
        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue
            
        # Determine tissue name from filename
        tissue_name = os.path.basename(file_path).split('.')[0].replace('gene_tpm_v10_', '')
        
        # Process file; pass the control marker gene of choice
        results_df, markers_df = process_gtex_file(
            file_path=file_path,
            target_gene_id=CPD_ID,
            gene_id_to_name={},
            control_marker_gene=CONTROL_MARKER
        )
        
        if results_df is not None and not results_df.empty:
            # Store correlation results
            all_results[tissue_name] = results_df
            
            # Plot correlation bar plot
            plot_angiogenesis_correlations(results_df, tissue_name)
            
            # Plot scatter plot with CPD vs. selected markers (PECAM1 and control marker)
            if markers_df is not None:
                plot_scatter_markers(markers_df, tissue_name)
            
            # Print significant correlations
            sig_df = results_df[results_df['fdr'] < 0.05]
            print(f"\nSignificant correlations in {tissue_name}:")
            print(sig_df[['gene_name', 'correlation', 'fdr']])
    
    # Optional: Save all correlation results to CSV
    if all_results:
        combined_df = pd.concat(
            {k: v.set_index('gene_name') for k, v in all_results.items()},
            axis=1
        )
        combined_df.to_csv("cpd_angiogenesis_correlations.csv")
