from gprofiler import GProfiler
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # Import seaborn
import numpy as np
import textwrap
from statsmodels.stats.multitest import multipletests
# Clean and prepare gene list

# Load the gene expression data
df = pd.read_csv("/Users/atharva/Downloads/DESeq2_results_KO_vs_WT.csv", sep=",")

#df = pd.read_csv("/Users/atharva/Downloads/rna_seq_with_(split eruption).txt", sep="\t")

# Clean column names
df.columns = df.columns.str.strip()
print(df.columns)
# Ensure correct column names
expected_columns = {"gene_name", "log2FC", "FDR"}
if not expected_columns.issubset(df.columns):
    raise ValueError(f"Missing required columns. Expected columns: {expected_columns}")

# Filter genes based on log2FC and FDR
pos_genes = df[(df["log2FC"] > 0.5) & (df["FDR"] < 0.05)]["gene_name"].dropna().unique().tolist()
print(len(pos_genes))
neg_genes = df[(df["log2FC"] < -0.5) & (df["FDR"] < 0.05)]["gene_name"].dropna().unique().tolist()
print(len(neg_genes))

# Initialize gProfiler
gp = GProfiler(return_dataframe=True)

# Run GO enrichment analysis
results = gp.profile(
    organism='mmusculus',  # Mouse genes (use 'hsapiens' for human)
    query=neg_genes,
    sources=['GO:BP'],     # GO Biological Processes
    significance_threshold_method='fdr',
    no_evidences=False,    # Include gene intersections
    user_threshold=0.05
)

# Process results
if not results.empty:
    # Compute FDR using Benjamini-Hochberg correction
    p_values = results['p_value'].values
    _, corrected_fdr, _, _ = multipletests(p_values, method='fdr_bh')

    # Add FDR to the results
    results['FDR'] = corrected_fdr

    # Filter significant results based on FDR
    significant_results = results[results['FDR'] < 0.05].copy()
    significant_results = results[results['significant']]

    # Compute -log10(p-value) for better visualization
    significant_results['neg_log_p'] = -np.log10(significant_results['p_value'])
    significant_results['neg_log_FDR'] = -np.log10(significant_results['FDR'])

    # Sort the results
    significant_results = significant_results.sort_values('FDR').head(12)
    significant_results['y_pos'] = range(len(significant_results) - 1, -1, -1)

    # Select relevant columns
    columns = ['source', 'native', 'name', 'p_value', 'FDR', 'precision', 
               'recall', 'intersection_size', 'intersections', 'neg_log_p','term_size', 'neg_log_FDR', 'y_pos','description']
    significant_results = significant_results[columns]
    # Remove GO terms linked to too many genes in the genome
    significant_results = significant_results[significant_results['term_size'] < 5000]

    # Print top 10 results
    print(f"Top {min(10, len(significant_results))} Significant GO Terms:")
    print(significant_results.head(10))

    significant_results = significant_results.sort_values('precision', ascending=False).head(12)
    significant_results['y_pos'] = range(len(significant_results) - 1, -1, -1)

    plt.figure(figsize=(14, 10))

    ax = plt.gca()

    # Scatter plot
    scatter = sns.scatterplot(
        data=significant_results,
        x="precision",
        y="y_pos",
        hue="neg_log_p",
        size="term_size",
        sizes=(25, 250),
        palette="viridis",
        edgecolor="black",
        legend="brief",
        ax=ax
    )

    # Improve text wrapping
    max_term_length = 90  # Characters per line for GO term names

    def wrap_text(text, max_length):
        return '\n'.join(textwrap.wrap(str(text), max_length))

    # Set y-axis labels with wrapped text
    wrapped_labels = [wrap_text(name, max_term_length) for name in significant_results["name"]]
    ax.set_yticks(significant_results["y_pos"])
    ax.set_yticklabels(wrapped_labels, fontsize=10)


    # Improve plot aesthetics
    ax.set_title("Top GO Terms ", fontsize=14, pad=20)
    ax.set_xlabel("Gene Ratio", fontsize=12)
    ax.set_ylabel("GO Term", fontsize=12)
    ax.grid(True, alpha=0.3)

    # Move legend
    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0.
    )

    plt.tight_layout()
    plt.show()

else:
    print("No significant GO terms found at FDR < 0.05")