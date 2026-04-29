import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))

# ============================================
# FILTER TO DOMINANT ISOFORMS ONLY
# ============================================

results_df = results_df[results_df["is_dominant"] == True]
print(f"After filtering to dominant isoforms: {len(results_df)} gene-sample pairs")
print(f"Unique genes: {results_df['gene_id'].nunique()}")

# ============================================
# CALCULATE SWITCHING VS STABLE
# ============================================

q2_df = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

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
# PLOT: Median TPM comparison — stable vs switching
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
ax1.set_title("Expression Level:\nStable vs Switching Genes (Dominant Isoforms Only)",
              fontsize=13)
ax1.yaxis.grid(True, linestyle="--", alpha=0.5)
ax1.set_axisbelow(True)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# --- Box plot of TPM distribution ---
ax2 = axes[1]

ax2.boxplot([stable["avg_total_tpm"], switching["avg_total_tpm"]],
            labels=[f"Stable genes\n(n={len(stable):,})",
                    f"Switching genes\n(n={len(switching):,})"],
            patch_artist=True,
            boxprops=dict(facecolor="#1D9E75", alpha=0.85),
            medianprops=dict(color="white", linewidth=2),
            flierprops=dict(marker="o", markersize=2, alpha=0.3))

boxes = ax2.patches
if len(boxes) > 1:
    boxes[1].set_facecolor("#E24B4A")

ax2.set_ylabel("log10(TPM + 1)", fontsize=12)
ax2.set_yscale("log")
ax2.set_title("Expression Distribution:\nStable vs Switching Genes (Dominant Isoforms Only)",
              fontsize=13)
ax2.yaxis.grid(True, linestyle="--", alpha=0.5)
ax2.set_axisbelow(True)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout()
output_path = os.path.expanduser(
    "~/finalproject/results/plots/plot_tpm_stable_vs_switching_dominant_only.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"\nPlot saved to {output_path}")
