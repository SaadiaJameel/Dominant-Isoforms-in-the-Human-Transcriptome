import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================
# CONFIGURATION
# ============================================
RESULTS_FILE = os.path.expanduser("~/finalproject/results/major_isoforms_tpm5.csv")
MIN_TPM = 5  # just for the plot title

print("Loading results...")
results_df = pd.read_csv(RESULTS_FILE)
print(f"Loaded {len(results_df)} gene-sample pairs")
print(f"Unique genes: {results_df['gene_id'].nunique()}")

# ============================================
# CALCULATE SWITCHING GENES
# ============================================

q2 = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

switching_genes = q2[q2["n_unique_major"] > 1].sort_values(
    "n_unique_major", ascending=False)
top20 = switching_genes.head(20)

print(f"\nTop 20 most switching genes:")
print(top20[["gene_name", "n_unique_major"]].to_string())

# ============================================
# PLOT
# ============================================

fig, ax = plt.subplots(figsize=(12, 8))

bars = ax.barh(top20["gene_name"], top20["n_unique_major"],
               color="#E24B4A", edgecolor="white", alpha=0.8)

for bar, val in zip(bars, top20["n_unique_major"]):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
            f"{val}", ha="left", va="center", fontsize=10)

ax.set_xlabel("Number of Unique Major Isoforms Across 18 Samples", fontsize=13)
ax.set_ylabel("Gene", fontsize=13)
ax.set_title(f"Top 20 Genes with Most Frequent Major Isoform Switching\n(TPM >= {MIN_TPM})",
             fontsize=14)
ax.invert_yaxis()

plt.tight_layout()

output_path = os.path.expanduser(
    f"~/finalproject/results/plots/plot_top20_switching_tpm{MIN_TPM}.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"\nPlot saved to {output_path}")
