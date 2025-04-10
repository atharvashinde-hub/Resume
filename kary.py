import mygene
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN
import numpy as np



# Initialize mygene
mg = mygene.MyGeneInfo()

# Load data
#df = pd.read_csv("/Users/atharva/Downloads/rna_seq_with_FDR(no thresh).txt", sep="\t")
#df = pd.read_csv("/Users/atharva/Downloads/rna_seq_with_(split eruption).txt", sep="\t")
#df = pd.read_csv("/Users/atharva/Downloads/rna_seq_with_FDR_FC.txt", sep="\t")
df = pd.read_csv("/Users/atharva/Downloads/DESeq2_results_3KO_vs_2WT.csv", sep=",")
# Clean column names
df.columns = df.columns.str.strip()

# Ensure correct column names
expected_columns = {"gene_name", "log2FC", "FDR"}
if not expected_columns.issubset(df.columns):
    raise ValueError(f"Missing required columns. Expected columns: {expected_columns}")

# Filter genes based on log2FC and FDR
pos_genes = df[(df["log2FC"] > 0.5) & (df["FDR"] < 0.05)][["gene_name", "log2FC"]]
neg_genes = df[(df["log2FC"] < -0.5) & (df["FDR"] < 0.05)][["gene_name", "log2FC"]]

# Combine positive and negative genes
filtered_genes = pd.concat([pos_genes, neg_genes])

# Query MyGeneInfo for genomic positions
gene_names = filtered_genes["gene_name"].tolist()
result = mg.querymany(gene_names, scopes="symbol", fields="genomic_pos", species="mouse")

# Process results into a DataFrame
# Process results into a DataFrame
gene_data = []
for gene in result:
    if "genomic_pos" in gene and gene["genomic_pos"]:
        genomic_positions = gene["genomic_pos"]
        
        # Ensure genomic_positions is a list
        if not isinstance(genomic_positions, list):
            genomic_positions = [genomic_positions]

        for genomic_pos in genomic_positions:  # Iterate through positions
            chromosome = genomic_pos.get("chr", "")
            start = genomic_pos.get("start", "")
            end = genomic_pos.get("end", "")
            
            # Retrieve log2FC value
            fc_values = filtered_genes.loc[filtered_genes["gene_name"] == gene["query"], "log2FC"].values
            if len(fc_values) > 0:
                fc = fc_values[0]
            else:
                fc = None  # Default if no match found
            
            gene_data.append([gene["query"], chromosome, start, end, fc])

df_genes = pd.DataFrame(gene_data, columns=["Gene", "Chromosome", "Start", "End", "log2FC"])


# Save to CSV
df_genes.to_csv("gene_locations.csv", index=False)
print("Gene location data saved to gene_locations.csv")
print(df_genes)

# ===============================
# Visualization of Gene Locations
# ===============================

# Convert Chromosome labels to numerical values for plotting
chromosomes = {ch: i for i, ch in enumerate(sorted(df_genes["Chromosome"].unique()))}
df_genes["Chromosome_num"] = df_genes["Chromosome"].map(chromosomes)

# Normalize log2FC for color mapping
norm = mcolors.Normalize(vmin=df_genes["log2FC"].min(), vmax=df_genes["log2FC"].max())
cmap = cm.get_cmap("coolwarm")  # Colormap from blue (negative) to red (positive)
df_genes["color"] = df_genes["log2FC"].apply(lambda fc: cmap(norm(fc)))

# Scale features for better clustering
scaler = StandardScaler()
scaled_features = scaler.fit_transform(df_genes[["Start", "Chromosome_num", "log2FC"]])

from sklearn.cluster import AgglomerativeClustering

agg_clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=0.2)  # You can adjust threshold
clusters = agg_clustering.fit_predict(scaled_features)

df_genes["Cluster"] = clusters

# Assign clusters
df_genes["Cluster"] = clusters
print(df_genes[["Gene", "Chromosome", "Start", "log2FC", "Cluster"]])
print(df_genes[df_genes["Cluster"] != -1][["Gene", "Chromosome", "Start", "log2FC", "Cluster"]])

# Plot gene positions on chromosomes
fig, ax = plt.subplots(figsize=(10, 6))

for _, row in df_genes.iterrows():
    ax.plot(
        [row["Start"], row["End"]],
        [row["Chromosome_num"], row["Chromosome_num"]],
        color=row["color"],
        alpha=0.7,
        linewidth=2
    )
    ax.text(
        row["Start"], row["Chromosome_num"],
        f"{row['Gene']} ({row['log2FC']:.2f})",
        fontsize=6, verticalalignment='bottom', alpha=0.7
    )

# Add color bar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Empty array for color bar
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("log2 Fold Change (log2FC)")

ax.set_yticks(list(chromosomes.values()))
ax.set_yticklabels(list(chromosomes.keys()))
ax.set_xlabel("Genomic Position (bp)")
ax.set_ylabel("Chromosome")
ax.set_title("Gene Locations on Chromosomes (Colored by Log2FC Spectrum)")
ax.set_title("Gene Locations on Chromosomes (Colored by Log2FC Spectrum)")
ax.set_title("COVERT anaysis (Chromosome Observation and Visuaization of Enriched RNA Transcripts)")

plt.show()

import matplotlib.patches as patches

# Plot setup
fig, ax = plt.subplots(figsize=(10, 6))

for _, row in df_genes.iterrows():
    ax.plot(
        [row["Start"], row["End"]],
        [row["Chromosome_num"], row["Chromosome_num"]],
        color=row["color"],
        alpha=0.7,
        linewidth=2
    )
    ax.text(
        row["Start"], row["Chromosome_num"],
        f"{row['Gene']} ({row['log2FC']:.2f})",
        fontsize=6, verticalalignment='bottom', alpha=0.7
    )

# Add color bar
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Empty array for color bar
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("log2 Fold Change (log2FC)")

# Loop through each chromosome
for chrom in df_genes["Chromosome_num"].unique():
    chrom_data = df_genes[df_genes["Chromosome_num"] == chrom]

    # Loop through each cluster within the chromosome
    for cluster_id in chrom_data["Cluster"].unique():
        if cluster_id == -1:  # Skip noise points
            continue

        # Get min and max Start positions for the cluster
        cluster_points = chrom_data[chrom_data["Cluster"] == cluster_id]
        min_start = cluster_points["Start"].min()
        max_start = cluster_points["Start"].max()
        chrom_pos = chrom  # Chromosome position on y-axis

        # Define rectangle
        rect = patches.Rectangle(
            (min_start, chrom_pos - 0.4),  # Bottom-left corner
            max_start - min_start,        # Width
            0.8,                          # Height
            linewidth=1.5,
            edgecolor="black",
            facecolor="none",
            linestyle="--"
        )
        ax.add_patch(rect)

# Labels and title
ax.set_xlabel("Genomic Position (bp)")
ax.set_ylabel("Chromosome")
ax.set_title("Gene Clusters with Bounding Boxes (Chromosome-Specific)")

plt.show()

# ... [Previous code remains unchanged until clustering section] ...

# ===============================
# Improved Clustering with Adjustable Parameters
# ===============================

# Set adjustable parameters
start_weight = 2.0    # Weight for genomic position (higher = more emphasis)
fc_weight = 1.0       # Weight for fold change
eps_value = 0.8       # DBSCAN neighborhood radius
min_samples_value = 2  # Minimum genes per cluster

# Initialize cluster column
df_genes["Cluster"] = -1  # -1 indicates noise

# Cluster within each chromosome separately
for chrom in df_genes["Chromosome"].unique():
    chrom_mask = df_genes["Chromosome"] == chrom
    chrom_data = df_genes[chrom_mask].copy()
    
    if len(chrom_data) < 2:  # Skip chromosomes with <2 genes
        continue
    
    # Extract and scale features
    X = chrom_data[["Start", "log2FC"]].values
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Apply feature weights
    X_weighted = X_scaled * np.array([start_weight, fc_weight])
    
    # Perform DBSCAN clustering
    dbscan = DBSCAN(eps=eps_value, min_samples=min_samples_value)
    clusters = dbscan.fit_predict(X_weighted)
    
    # Assign cluster IDs (add chromosome prefix for uniqueness)
    df_genes.loc[chrom_mask, "Cluster"] = [f"{chrom}_{cid}" if cid != -1 else -1 for cid in clusters]

# Filter meaningful clusters (≥2 genes)
cluster_counts = df_genes[df_genes["Cluster"] != -1].groupby("Cluster").size()
valid_clusters = cluster_counts[cluster_counts >= 2].index
df_genes["Cluster"] = df_genes["Cluster"].where(df_genes["Cluster"].isin(valid_clusters), -1)

# ===============================
# Enhanced Visualization
# ===============================

fig, ax = plt.subplots(figsize=(12, 8))

# Plot genes
for _, row in df_genes.iterrows():
    ax.plot(
        [row["Start"], row["End"]],
        [row["Chromosome_num"], row["Chromosome_num"]],
        color=row["color"],
        alpha=0.7,
        linewidth=2
    )
    ax.text(
        row["Start"], row["Chromosome_num"],
        f"{row['Gene']} ({row['log2FC']:.2f})",
        fontsize=6, verticalalignment='bottom', alpha=0.7
    )

# Add cluster bounding boxes
for cluster_id in valid_clusters:
    chrom = cluster_id.split("_")[0]
    cluster_data = df_genes[(df_genes["Cluster"] == cluster_id) & (df_genes["Chromosome"] == chrom)]
    
    if len(cluster_data) < 2:
        continue
    
    min_start = cluster_data["Start"].min()
    max_start = cluster_data["Start"].max()
    chrom_pos = chromosomes[chrom]
    
    rect = patches.Rectangle(
        (min_start, chrom_pos - 0.4),
        max_start - min_start,
        0.8,
        linewidth=1.5,
        edgecolor="black",
        facecolor="none",
        linestyle="--"
    )
    ax.add_patch(rect)

# Add color bar and labels
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label("log2 Fold Change (log2FC)")

ax.set_yticks(list(chromosomes.values()))
ax.set_yticklabels(list(chromosomes.keys()))
ax.set_xlabel("Genomic Position (bp)")
ax.set_ylabel("Chromosome")
ax.set_title("COVERT Analysis with Improved Clustering")

plt.show()