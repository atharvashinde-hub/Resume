import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
from pyvis.network import Network
import requests  # For API calls
import time  # To handle API rate limits

# --------------------------- Step 1: Load and Preprocess Data --------------------------- #
df = pd.read_csv("/Users/atharva/Downloads/rna_seq_with_FDR(no thresh).txt", sep="\t")
df.columns = df.columns.str.strip()

expected_columns = {"gene_name", "log2FC", "FDR"}
if not expected_columns.issubset(df.columns):
    raise ValueError(f"Missing columns! Expected: {expected_columns}, Found: {df.columns}")

df = df.dropna(subset=["gene_name", "log2FC", "FDR"])

# --------------------------- Step 2: Filter DEGs --------------------------- #
# Filter DEGs with FDR < 0.05
filtered_deg = df[df['FDR'] < 0.1].copy()

# --------------------------- Step 3: Load TF Data from Multiple Sources --------------------------- #
# Load new TF-target data and standardize names to uppercase
new_tf_targets = []
with open("/Users/atharva/Downloads/fibroblast.TFs.gmt", "r") as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) < 2:
            continue  # Skip lines with no targets
        tf = parts[0].upper()  # Convert TF to uppercase
        targets = [t.upper() for t in parts[1:]]  # Convert targets to uppercase
        for target in targets:
            new_tf_targets.append({'TF': tf, 'Target': target, 'Interaction': 'Unknown', 'PMID': 'None', 'Source': 'new_tf_targets'})

new_tf_targets = pd.DataFrame(new_tf_targets)

# Load mouse.source
mouse_source = pd.read_csv("/Users/atharva/Downloads/mouse.source", sep="\t", header=None,
                           names=['TF', 'ID', 'Target', 'PMID'])
mouse_source = mouse_source[['TF', 'Target', 'PMID']]
mouse_source['Interaction'] = 'Unknown'  # No interaction type provided
mouse_source['Source'] = 'mouse.source'

# Load trrust_rawdata.mouse.tsv
trrust_data = pd.read_csv("/Users/atharva/Downloads/trrust_rawdata.mouse.tsv", sep="\t", header=None,
                          names=['TF', 'Target', 'Interaction', 'PMID'])
trrust_data['Source'] = 'trrust'

# Standardize gene names in all TF datasets
for dataset in [new_tf_targets, mouse_source, trrust_data]:
    dataset['TF'] = dataset['TF'].str.upper()
    dataset['Target'] = dataset['Target'].str.upper()

# Combine all TF datasets
tf_df = pd.concat([new_tf_targets, mouse_source, trrust_data], ignore_index=True)

# Standardize gene names in DEG list and expression data
filtered_deg['gene_name'] = filtered_deg['gene_name'].str.upper()
df_aggregated = df.groupby('gene_name').agg({'log2FC': 'mean', 'FDR': 'mean'}).reset_index()
df_aggregated['gene_name'] = df_aggregated['gene_name'].str.upper()

# Check TF overlap
all_tfs = tf_df['TF'].unique()
overlap_tfs = set(all_tfs).intersection(set(filtered_deg['gene_name']))
print(f"Number of overlapping TFs: {len(overlap_tfs)}")
print("Overlapping TFs:", overlap_tfs)

# --------------------------- Step 4: Identify TF DEGs --------------------------- #
tf_deg = filtered_deg[filtered_deg['gene_name'].isin(all_tfs)]

# --------------------------- Step 5: Find Targets & Their Log2FC --------------------------- #
gene_expression = df_aggregated.set_index('gene_name')[['log2FC', 'FDR']].to_dict('index')

result = []
for _, row in tf_deg.iterrows():
    tf = row['gene_name']
    targets = tf_df[tf_df['TF'] == tf]
    for _, target_row in targets.iterrows():
        target_gene = target_row['Target']
        target_info = gene_expression.get(target_gene, {})
        if target_info.get('FDR', 1) < 0.05 and abs(target_info.get('log2FC', 0)) >= 0.5:
            result.append({
                'TF': tf,
                'TF_log2FC': row['log2FC'],
                'TF_FDR': row['FDR'],
                'Target': target_gene,
                'Interaction_Type': target_row['Interaction'],
                'PMID': target_row['PMID'],
                'Target_log2FC': target_info.get('log2FC', np.nan),
                'Target_FDR': target_info.get('FDR', np.nan),
                'Source': target_row['Source']
            })

# --------------------------- Step 6: Create Final DataFrame --------------------------- #
final_df = pd.DataFrame(result)

if not final_df.empty:
    final_df = final_df.dropna(subset=['Target_log2FC'])
    # --------------------------- Visualization --------------------------- #
    plt.figure(figsize=(12, 8))
    sns.scatterplot(data=final_df, x='TF_log2FC', y='Target_log2FC', hue='Interaction_Type', style='Source', s=100, alpha=0.7)
    plt.axhline(0, color='grey', linestyle='--')
    plt.axvline(0, color='grey', linestyle='--')
    plt.title('TF-Target Expression Relationships')
    plt.show()

    # --------------------------- Save Results --------------------------- #
    final_df.to_csv('/Users/atharva/Downloads/tf_target_analysis_combined.csv', index=False, sep="\t")
    print(f"Found {len(tf_deg)} TFs among DEGs")
    print(f"Identified {len(final_df)} TF-target relationships with target FDR < 0.05 and |log2FC| >= 0.5")
else:
    print("No TF-target relationships found. Exiting.")

# --------------------------- Network Visualization --------------------------- #
G = nx.DiGraph()

# Add nodes (TFs and Targets)
for _, row in final_df.iterrows():
    tf = row['TF']
    target = row['Target']
    
    # Add TF and target nodes
    if tf not in G:
        G.add_node(tf, log2FC=row['TF_log2FC'])
    if target not in G:
        G.add_node(target, log2FC=row['Target_log2FC'])
    
    # Determine edge color based on FC (activation or inhibition)
    if row['TF_log2FC'] > 0:
        edge_color = "red"  # Activation
    else:
        edge_color = "blue"  # Inhibition

    # Add edge
    G.add_edge(tf, target, color=edge_color)

# Interactive Visualization
net = Network(notebook=True, directed=True)

# Function to fetch gene description from Ensembl REST API
def get_gene_description(gene_name):
    try:
        # Ensembl REST API endpoint
        url = f"https://rest.ensembl.org/lookup/symbol/mus_musculus/{gene_name}?expand=1"
        headers = {"Content-Type": "application/json"}
        response = requests.get(url, headers=headers)
        
        if response.ok:
            data = response.json()
            return data.get('description', 'No description available')
        else:
            return 'Description not found'
    except Exception as e:
        return f"API error: {str(e)}"

# Add nodes with gene descriptions
for node, data in G.nodes(data=True):
    log2FC = data.get("log2FC", 0)
    color = "red" if log2FC > 0 else "blue"
    size = abs(log2FC) * 10 + 10
    
    # Fetch gene description
    description = get_gene_description(node)
    
    # Create hover text
    hover_text = f"Gene: {node}<br>log2FC: {log2FC:.2f}<br>Description: {description}"
    
    # Add node to network
    net.add_node(node, title=hover_text, color=color, size=size)

    # Add a delay to avoid hitting API rate limits
    time.sleep(0.1)

# Add edges
for source, target, data in G.edges(data=True):
    net.add_edge(source, target, color=data["color"])

# Save and display
net.show("tf_target_network.html")