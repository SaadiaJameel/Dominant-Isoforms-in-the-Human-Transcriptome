import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
q2_df = pd.read_csv(os.path.expanduser("~/finalproject/results/q2_isoform_consistency.csv"))

# Get MANE agreement per gene averaged across all samples
mane_per_gene = results_df.dropna(subset=["matches_mane"]).groupby("gene_id").agg(
    mane_agreement=("matches_mane", "mean")
).reset_index()
mane_per_gene["mane_agreement_pct"] = mane_per_gene["mane_agreement"] * 100

# Merge with switching info
merged = q2_df.merge(mane_per_gene, on="gene_id", how="inner")

# Split into two groups
stable = merged[merged["n_unique_major"] == 1]
switching = merged[merged["n_unique_major"] > 1]

stable_mane = stable["mane_agreement_pct"].mean()
switching_mane = switching["mane_agreement_pct"].mean()

print(f"Stable genes:    n={len(stable)},    MANE agreement = {stable_mane:.1f}%")
print(f"Switching genes: n={len(switching)}, MANE agreement = {switching_mane:.1f}%")

# ── Plot ──────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 6))

categories = [f"Stable genes\n(n={len(stable):,})", 
              f"Switching genes\n(n={len(switching):,})"]
values = [stable_mane, switching_mane]
colors = ["#1D9E75", "#E24B4A"]

bars = ax.bar(categories, values, color=colors, alpha=0.85,
              edgecolor="white", width=0.5)

# Value labels on bars
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
            f"{val:.1f}%", ha="center", va="bottom",
            fontsize=14, fontweight="bold")

# Overall average line
overall = results_df.dropna(subset=["matches_mane"])["matches_mane"].mean() * 100
ax.axhline(y=overall, color="black", linestyle="--", linewidth=1.5,
           label=f"Overall average ({overall:.1f}%)")

ax.set_ylabel("MANE_Select agreement rate (%)", fontsize=13)
ax.set_title("MANE_Select Agreement:\nStable vs Switching Genes", fontsize=14)
ax.set_ylim(0, 100)
ax.legend(fontsize=11)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
output_path = os.path.expanduser(
    "~/finalproject/results/plots/plot_mane_stable_vs_switching.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"Plot saved to {output_path}")
