import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
print(f"Loaded {len(results_df)} rows")

def categorize(pct):
    if pct == 100:
        return "100% dominant"
    elif pct > 50:
        return "50-99% dominant"
    else:
        return "< 50% dominant"

results_df["dominance_cat"] = results_df["major_pct"].apply(categorize)

counts = results_df.groupby(["sample", "dominance_cat"]).size().unstack(fill_value=0)

for col in ["< 50% dominant", "50-99% dominant", "100% dominant"]:
    if col not in counts.columns:
        counts[col] = 0
counts = counts[["< 50% dominant", "50-99% dominant", "100% dominant"]]
counts = counts.sort_index()

fig, ax = plt.subplots(figsize=(18, 7))

samples = counts.index.tolist()
x = np.arange(len(samples))
width = 0.25

colors = ["#E24B4A", "#378ADD", "#1D9E75"]
labels = ["< 50% (no dominant isoform)", "50-99% (dominant isoform)", "100% (single isoform only)"]

ax.bar(x - width, counts["< 50% dominant"],  width, label=labels[0], color=colors[0], alpha=0.85)
ax.bar(x,          counts["50-99% dominant"], width, label=labels[1], color=colors[1], alpha=0.85)
ax.bar(x + width,  counts["100% dominant"],   width, label=labels[2], color=colors[2], alpha=0.85)

ax.set_xlabel("Sample", fontsize=13)
ax.set_ylabel("Number of protein-coding genes", fontsize=13)
ax.set_title("Isoform Dominance Categories Across 18 Geuvadis Samples", fontsize=15, pad=15)
ax.set_xticks(x)
ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
ax.legend(fontsize=11, loc="upper right")
ax.set_ylim(0, counts.max().max() * 1.2)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
output_path = os.path.expanduser("~/finalproject/results/plots/plot_dominance_3bars.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"Plot saved to {output_path}")
