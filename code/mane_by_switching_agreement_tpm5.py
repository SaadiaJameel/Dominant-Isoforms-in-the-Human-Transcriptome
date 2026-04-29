import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms_tpm5.csv"))
q2_df = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

print(f"Loaded {len(results_df)} rows")
print(f"Loaded {len(q2_df)} genes from Q2 analysis")

# ============================================
# MERGE SWITCHING INFO WITH MANE AGREEMENT
# ============================================

# q2_df has gene_id and n_unique_major (number of unique major isoforms)
# results_df has matches_mane per gene per sample

# Get MANE agreement per gene (average across all samples)
mane_per_gene = results_df.dropna(subset=["matches_mane"]).groupby("gene_id").agg(
    mane_agreement=("matches_mane", "mean"),
    n_samples=("sample", "count")
).reset_index()
mane_per_gene["mane_agreement_pct"] = mane_per_gene["mane_agreement"] * 100

# Merge with switching info
merged = q2_df.merge(mane_per_gene, on="gene_id", how="inner")
print(f"\nGenes with both switching info and MANE data: {len(merged)}")

# ============================================
# CALCULATE MANE AGREEMENT BY N_UNIQUE_MAJOR
# ============================================

# Cap at 10+ for display
def cap_switching(n):
    if n >= 10:
        return "10+"
    return str(int(n))

merged["switching_group"] = merged["n_unique_major"].apply(cap_switching)

summary = merged.groupby("switching_group").agg(
    n_genes=("mane_agreement_pct", "count"),
    mane_agreement=("mane_agreement_pct", "mean")
).reset_index()

# Sort properly
order = [str(i) for i in range(1, 10)] + ["10+"]
summary["switching_group"] = pd.Categorical(
    summary["switching_group"], categories=order, ordered=True)
summary = summary.sort_values("switching_group")

print("\nMANE agreement by switching category:")
print(summary[["switching_group", "n_genes", "mane_agreement"]].to_string(index=False))

# ============================================
# PLOT
# ============================================

fig, ax = plt.subplots(figsize=(13, 6))

x = np.arange(len(summary))

# Color bars by switching level
colors = ["#1D9E75" if i == 0 else
          "#378ADD" if i <= 2 else
          "#EF9F27" if i <= 5 else
          "#E24B4A"
          for i in range(len(summary))]

bars = ax.bar(x, summary["mane_agreement"],
              color=colors, alpha=0.85, edgecolor="white")

# Value labels
for bar, val in zip(bars, summary["mane_agreement"]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{val:.1f}%", ha="center", va="bottom", fontsize=9.5)

# Overall MANE agreement line
overall = results_df.dropna(subset=["matches_mane"])["matches_mane"].mean() * 100
ax.axhline(y=overall, color="black", linestyle="--", linewidth=1.5,
           label=f"Overall MANE agreement ({overall:.1f}%)")

ax.set_xlabel("Number of unique major isoforms across 18 samples", fontsize=12)
ax.set_ylabel("MANE_Select agreement rate (%)", fontsize=12)
ax.set_title("MANE_Select Agreement Rate by Isoform Switching Level",
             fontsize=14, pad=15)
ax.set_xticks(x)
ax.set_xticklabels(summary["switching_group"], fontsize=11)
ax.set_ylim(0, 100)
ax.legend(fontsize=11)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Color legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#1D9E75", alpha=0.85, label="Perfectly stable (x=1)"),
    Patch(facecolor="#378ADD", alpha=0.85, label="Minor switching (x=2-3)"),
    Patch(facecolor="#EF9F27", alpha=0.85, label="Moderate switching (x=4-6)"),
    Patch(facecolor="#E24B4A", alpha=0.85, label="Frequent switching (x=7+)"),
]
ax.legend(handles=legend_elements + [
    plt.Line2D([0], [0], color="black", linestyle="--",
               label=f"Overall MANE agreement ({overall:.1f}%)")
], fontsize=10, loc="upper right")

# Gene counts below x axis
for i, (_, row) in enumerate(summary.iterrows()):
    ax.text(i, -5, f"n={int(row['n_genes'])}",
            ha="center", va="top", fontsize=8, color="gray")

plt.tight_layout()
output_path = os.path.expanduser(
    "~/finalproject/results/plots/plot_mane_by_switching.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to {output_path}")
