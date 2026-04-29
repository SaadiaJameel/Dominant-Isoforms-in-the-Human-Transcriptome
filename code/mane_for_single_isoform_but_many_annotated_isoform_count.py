import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
always_100_genes = pd.read_csv(os.path.expanduser("~/finalproject/results/always_100_genes.csv"))

print(f"Loaded {len(results_df)} rows")
print(f"Loaded {len(always_100_genes)} always-100% genes")

# ============================================
# COPY ALL GENES
# ============================================

blue_genes = always_100_genes.copy()
print(f"\nTotal always-100% genes: {len(blue_genes)}")

# Drop genes with no MANE_Select
blue_genes_mane = blue_genes.dropna(subset=["mane_transcript"])
print(f"Of these, genes with a MANE_Select: {len(blue_genes_mane)}")

# ============================================
# CALCULATE MANE AGREEMENT RATE PER ISOFORM COUNT
# ============================================

def cap_isoforms(n):
    if n >= 15:
        return "15+"
    return str(int(n))

blue_genes_mane["isoform_group"] = blue_genes_mane["n_isoforms"].apply(cap_isoforms)

summary = blue_genes_mane.groupby("isoform_group").agg(
    n_genes=("matches_mane", "count"),
    mane_agreement=("matches_mane", "mean")
).reset_index()

summary["mane_agreement_pct"] = summary["mane_agreement"] * 100

order = [str(i) for i in range(1, 15)] + ["15+"]
summary["isoform_group"] = pd.Categorical(summary["isoform_group"], categories=order, ordered=True)
summary = summary.sort_values("isoform_group")

print("\nMANE agreement by isoform count:")
print(summary[["isoform_group", "n_genes", "mane_agreement_pct"]].to_string(index=False))

# ============================================
# OVERALL RATE FROM SAME DATA
# ============================================

overall = blue_genes_mane["matches_mane"].mean() * 100
print(f"\nOverall MANE agreement (always-100% genes): {overall:.1f}%")

# ============================================
# PLOT
# ============================================

fig, ax1 = plt.subplots(figsize=(13, 6))

x = np.arange(len(summary))

bars = ax1.bar(x, summary["mane_agreement_pct"],
               color="#378ADD", alpha=0.85, edgecolor="white", label="MANE agreement %")

for bar, val in zip(bars, summary["mane_agreement_pct"]):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f"{val:.1f}%", ha="center", va="bottom", fontsize=9)

ax1.axhline(y=overall, color="red", linestyle="--", linewidth=1.5,
            label=f"Overall MANE agreement ({overall:.1f}%)")

ax1.set_ylabel("MANE_Select agreement rate (%)", fontsize=12)
ax1.set_ylim(0, 100)
ax1.set_xticks(x)
ax1.set_xticklabels(summary["isoform_group"], fontsize=11)
ax1.set_xlabel("Number of annotated isoforms in GENCODE v37", fontsize=12)
ax1.set_title("MANE_Select Agreement Rate for Always-100% Dominant Genes\n"
              "by Number of Annotated Isoforms",
              fontsize=13)
ax1.legend(fontsize=11, loc="upper right")
ax1.yaxis.grid(True, linestyle="--", alpha=0.5)
ax1.set_axisbelow(True)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

for i, (_, row) in enumerate(summary.iterrows()):
    ax1.text(i, -6, f"n={int(row['n_genes'])}",
             ha="center", va="top", fontsize=8, color="gray")

plt.tight_layout()
output_path = os.path.expanduser("~/finalproject/results/plots/plot_mane_by_isoform_count.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to {output_path}")
