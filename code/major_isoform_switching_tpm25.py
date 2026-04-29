import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================
# CONFIGURATION
# ============================================
RESULTS_FILE = os.path.expanduser("~/finalproject/results/major_isoforms_tpm25.csv")
MIN_TPM = 25  # just for the plot title

print("Loading results...")
results_df = pd.read_csv(RESULTS_FILE)
print(f"Loaded {len(results_df)} gene-sample pairs")
print(f"Unique genes: {results_df['gene_id'].nunique()}")

# ============================================
# CALCULATE SWITCHING
# ============================================

q2 = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

# ============================================
# PLOT
# ============================================

fig, ax = plt.subplots(figsize=(10, 6))

unique_counts = q2["n_unique_major"].value_counts().sort_index()

bars = ax.bar(unique_counts.index, unique_counts.values,
              color="steelblue", edgecolor="white", alpha=0.8)

for bar, val in zip(bars, unique_counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
            f"{val:,}", ha="center", va="bottom", fontsize=9)

ax.set_xlabel("Number of Unique Major Isoforms Across 18 Samples", fontsize=13)
ax.set_ylabel("Number of Genes", fontsize=13)
ax.set_title(f"Distribution of Major Isoform Switching\nAcross Protein-Coding Genes (TPM >= {MIN_TPM})",
             fontsize=14)
ax.set_xticks(unique_counts.index)

plt.tight_layout()

output_path = os.path.expanduser(
    f"~/finalproject/results/plots/plot_switching_distribution_tpm{MIN_TPM}.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"\nPlot saved to {output_path}")
