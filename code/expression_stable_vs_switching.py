import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
q2_df = pd.read_csv(os.path.expanduser("~/finalproject/results/q2_isoform_consistency.csv"))

# ============================================
# GET AVERAGE TPM PER GENE ACROSS ALL SAMPLES
# ============================================

avg_tpm = results_df.groupby(["gene_id", "gene_name"]).agg(
    avg_total_tpm=("total_tpm", "mean")
).reset_index()

# Merge with switching info
merged = q2_df.merge(avg_tpm, on="gene_id", how="inner")

# Split into stable vs switching
stable = merged[merged["n_unique_major"] == 1]
switching = merged[merged["n_unique_major"] > 1]

print(f"Stable genes:    n={len(stable)},   median TPM = {stable['avg_total_tpm'].median():.2f}")
print(f"Switching genes: n={len(switching)}, median TPM = {switching['avg_total_tpm'].median():.2f}")
print(f"Mean TPM stable:    {stable['avg_total_tpm'].mean():.2f}")
print(f"Mean TPM switching: {switching['avg_total_tpm'].mean():.2f}")

# ============================================
# PLOT 1: Median TPM comparison — stable vs switching
# ============================================

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# --- Bar chart of median TPM ---
ax1 = axes[0]

categories = [f"Stable genes\n(n={len(stable):,})",
              f"Switching genes\n(n={len(switching):,})"]
medians = [stable["avg_total_tpm"].median(),
           switching["avg_total_tpm"].median()]
colors = ["#1D9E75", "#E24B4A"]

bars = ax1.bar(categories, medians, color=colors, alpha=0.85,
               edgecolor="white", width=0.5)

for bar, val in zip(bars, medians):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f"{val:.1f}", ha="center", va="bottom",
             fontsize=13, fontweight="bold")

ax1.set_ylabel("Median total gene TPM across 18 samples", fontsize=12)
ax1.set_title("Expression Level:\nStable vs Switching Genes", fontsize=13)
ax1.yaxis.grid(True, linestyle="--", alpha=0.5)
ax1.set_axisbelow(True)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# --- Box plot for distribution ---
ax2 = axes[1]

# Log scale for better visualization
stable_log = np.log10(stable["avg_total_tpm"] + 1)
switching_log = np.log10(switching["avg_total_tpm"] + 1)

bp = ax2.boxplot([stable_log, switching_log],
                 patch_artist=True,
                 medianprops=dict(color="white", linewidth=2),
                 whiskerprops=dict(linewidth=1.5),
                 capprops=dict(linewidth=1.5),
                 flierprops=dict(marker="o", markersize=2, alpha=0.3))

bp["boxes"][0].set_facecolor("#1D9E75")
bp["boxes"][0].set_alpha(0.85)
bp["boxes"][1].set_facecolor("#E24B4A")
bp["boxes"][1].set_alpha(0.85)

ax2.set_xticks([1, 2])
ax2.set_xticklabels([f"Stable genes\n(n={len(stable):,})",
                     f"Switching genes\n(n={len(switching):,})"], fontsize=11)
ax2.set_ylabel("log10(TPM + 1)", fontsize=12)
ax2.set_title("Expression Distribution:\nStable vs Switching Genes", fontsize=13)
ax2.yaxis.grid(True, linestyle="--", alpha=0.5)
ax2.set_axisbelow(True)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.suptitle("Do Switching Genes Have Lower Expression in LCLs?",
             fontsize=14, fontweight="bold", y=1.02)

plt.tight_layout()
output_path = os.path.expanduser(
    "~/finalproject/results/plots/plot_expression_stable_vs_switching.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to {output_path}")

# ============================================
# ALSO SHOW TOP 20 STABLE GENES WITH HIGHEST TPM
# ============================================

print("\n--- Top 20 most stable AND highly expressed genes ---")
top_stable = stable.sort_values("avg_total_tpm", ascending=False).head(20)
print(top_stable[["gene_id", "avg_total_tpm"]].to_string(index=False))

print("\n--- Top 20 most switching genes with their TPM ---")
top_switching = switching.sort_values("n_unique_major", ascending=False).head(20)
print(top_switching[["gene_id", "n_unique_major", "avg_total_tpm"]].to_string(index=False))
