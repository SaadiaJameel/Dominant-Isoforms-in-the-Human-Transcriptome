import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
print(f"Loaded {len(results_df)} rows")

# Get genes that are 100% dominant in each sample
genes_100 = results_df[results_df["major_pct"] == 100]

# For each gene, get the list of samples where it is 100% dominant
gene_sample_membership = genes_100.groupby("gene_name")["sample"].apply(list).reset_index()

# Build the membership list for upsetplot
memberships = gene_sample_membership["sample"].tolist()

# Create the upset data
upset_data = from_memberships(memberships)

# Plot
fig = plt.figure(figsize=(20, 8))
upset = UpSet(
    upset_data,
    subset_size="count",
    show_counts=True,
    sort_by="cardinality",
    min_subset_size=50
)
upset.plot(fig)

plt.suptitle("Overlap of 100% Dominant Genes Across 18 Geuvadis Samples",
             fontsize=14, y=1.02)

output_path = os.path.expanduser("~/finalproject/results/plots/upset_100pct_genes.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot saved to {output_path}")

# Also print a quick summary
n_all_18 = sum(1 for m in memberships if len(m) == 18)
n_at_least_15 = sum(1 for m in memberships if len(m) >= 15)
print(f"\nGenes 100% dominant in all 18 samples:     {n_all_18}")
print(f"Genes 100% dominant in at least 15 samples: {n_at_least_15}")
print(f"Total genes ever 100% dominant:             {len(memberships)}")
