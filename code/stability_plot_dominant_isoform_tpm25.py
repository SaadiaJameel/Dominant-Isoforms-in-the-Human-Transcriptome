import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================
# CONFIGURATION — change to compare versions
# ============================================
RESULTS_FILE = os.path.expanduser("~/finalproject/results/major_isoforms_tpm25.csv")
MIN_TPM = 25 # just for the plot title

print("Loading results...")
results_df = pd.read_csv(RESULTS_FILE)
print(f"Loaded {len(results_df)} gene-sample pairs")
print(f"Unique genes: {results_df['gene_id'].nunique()}")

# ============================================
# FILTER TO DOMINANT ISOFORMS ONLY
# ============================================

results_df = results_df[results_df["is_dominant"] == True]
print(f"\nAfter filtering to dominant isoforms: {len(results_df)} gene-sample pairs")
print(f"Unique genes: {results_df['gene_id'].nunique()}")

# ============================================
# CALCULATE SWITCHING VS STABLE
# ============================================

q2 = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

n_consistent = (q2["n_unique_major"] == 1).sum()
n_switching  = (q2["n_unique_major"] > 1).sum()
n_total      = len(q2)

pct_consistent = n_consistent / n_total * 100
pct_switching  = n_switching  / n_total * 100

print(f"\nTotal expressed genes: {n_total:,}")
print(f"Consistent (stable):   {n_consistent:,} ({pct_consistent:.1f}%)")
print(f"Switching:             {n_switching:,}  ({pct_switching:.1f}%)")

# ============================================
# PLOT
# ============================================

fig, ax = plt.subplots(figsize=(8, 8))

wedges, texts, autotexts = ax.pie(
    [n_consistent, n_switching],
    labels=[f"Consistent\n({n_consistent:,} genes)",
            f"Switching\n({n_switching:,} genes)"],
    autopct="%1.1f%%",
    colors=["#378ADD", "#E24B4A"],
    startangle=90,
    textprops={"fontsize": 12}
)

for autotext in autotexts:
    autotext.set_fontsize(13)
    autotext.set_fontweight("bold")

ax.set_title(f"Stability of Dominant Isoform Identity\nAcross All 18 Samples (TPM >= {MIN_TPM})",
             fontsize=14, pad=20)

plt.tight_layout()

output_path = os.path.expanduser(
    f"~/finalproject/results/plots/plot_dominant_stability_pie_tpm{MIN_TPM}.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"\nPlot saved to {output_path}")
