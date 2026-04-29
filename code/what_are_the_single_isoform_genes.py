import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
print(f"Loaded {len(results_df)} rows")

# ============================================
# GET THE GENES ALWAYS 100% DOMINANT
# ============================================

always_100 = results_df.groupby("gene_id").filter(
    lambda x: (x["major_pct"] == 100).all()
)
always_100_genes = always_100[["gene_id", "gene_name",
                                "major_transcript",
                                "mane_transcript",
                                "matches_mane"]].drop_duplicates(subset="gene_id")

print(f"\nGenes always 100% dominant across all 18 samples: {len(always_100_genes)}")

# ============================================
# QUESTION 1: MANE_Select agreement rate
# ============================================

print("\n--- Q1: MANE_Select agreement ---")
has_mane = always_100_genes.dropna(subset=["mane_transcript"])
mane_agreement = has_mane["matches_mane"].mean() * 100
print(f"Genes with a MANE_Select transcript: {len(has_mane)}")
print(f"MANE_Select agreement rate:          {mane_agreement:.1f}%")

overall_mane = results_df.dropna(subset=["matches_mane"])
overall_rate = overall_mane["matches_mane"].mean() * 100
print(f"Overall MANE agreement rate:         {overall_rate:.1f}%")
print(f"Difference:                          +{mane_agreement - overall_rate:.1f}%")

# ============================================
# QUESTION 2: Number of annotated isoforms
# ============================================

print("\n--- Q2: Number of annotated isoforms ---")

gtf_file = "/opt/ccb/data/annotations/gencode.v37.annotation.gtf"
transcript_to_gene = {}
print("Parsing GTF for isoform counts (this may take a minute)...")

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "transcript":
            continue
        attributes = fields[8]
        gene_id = ""
        transcript_id = ""
        gene_type = ""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1].split(".")[0]
            if attr.startswith("transcript_id"):
                transcript_id = attr.split('"')[1].split(".")[0]
            if attr.startswith("gene_type"):
                gene_type = attr.split('"')[1]
        if gene_type == "protein_coding" and gene_id and transcript_id:
            if gene_id not in transcript_to_gene:
                transcript_to_gene[gene_id] = []
            transcript_to_gene[gene_id].append(transcript_id)

isoform_counts = {gene: len(transcripts)
                  for gene, transcripts in transcript_to_gene.items()}
isoform_df = pd.DataFrame(list(isoform_counts.items()),
                           columns=["gene_id", "n_isoforms"])

always_100_with_counts = always_100_genes.merge(isoform_df, on="gene_id", how="left")

median_always = always_100_with_counts["n_isoforms"].median()
median_all    = isoform_df["n_isoforms"].median()
mean_always   = always_100_with_counts["n_isoforms"].mean()
mean_all      = isoform_df["n_isoforms"].mean()

print(f"Median isoforms for always-100% genes:       {median_always:.1f}")
print(f"Median isoforms for all protein-coding genes: {median_all:.1f}")
print(f"Mean isoforms for always-100% genes:         {mean_always:.1f}")
print(f"Mean isoforms for all protein-coding genes:  {mean_all:.1f}")

print("\nIsoform count distribution for always-100% genes:")
print(always_100_with_counts["n_isoforms"].value_counts().sort_index().head(10).to_string())

# ============================================
# QUESTION 3: First 20 gene names
# ============================================

print("\n--- Q3: Sample of always-100% genes ---")
print(always_100_genes["gene_name"].head(20).to_string())

# ============================================
# PLOT 1: Mean and Median comparison bar chart
# ============================================

print("\nGenerating Plot 1: Mean and Median comparison...")

fig, ax = plt.subplots(figsize=(8, 6))

categories = ["Always 100%\ndominant genes", "All protein-coding\ngenes"]
medians = [median_always, median_all]
means   = [mean_always, mean_all]

x = np.arange(len(categories))
width = 0.35

bars1 = ax.bar(x - width/2, medians, width,
               label="Median", color="#378ADD", alpha=0.85)
bars2 = ax.bar(x + width/2, means, width,
               label="Mean", color="#1D9E75", alpha=0.85)

# Add value labels on bars
for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            f"{bar.get_height():.1f}", ha="center", va="bottom",
            fontsize=11, fontweight="bold")
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
            f"{bar.get_height():.1f}", ha="center", va="bottom",
            fontsize=11, fontweight="bold")

ax.set_ylabel("Number of annotated isoforms", fontsize=13)
ax.set_title("Annotated Isoform Count:\nAlways-100% Genes vs All Protein-Coding Genes",
             fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=12)
ax.legend(fontsize=11)
ax.set_ylim(0, max(means) * 1.25)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
path1 = os.path.expanduser("~/finalproject/results/plots/plot_isoform_count_comparison.png")
plt.savefig(path1, dpi=150)
plt.close()
print(f"  Plot 1 saved to {path1}")

# ============================================
# PLOT 2: Isoform count distribution for always-100% genes
# ============================================

print("Generating Plot 2: Isoform count distribution...")

dist = always_100_with_counts["n_isoforms"].value_counts().sort_index()

# Only show up to 15 isoforms for clarity, group the rest as "15+"
dist_trimmed = dist[dist.index <= 15].copy()
rest = dist[dist.index > 15].sum()
if rest > 0:
    dist_trimmed[">15"] = rest

fig, ax = plt.subplots(figsize=(12, 6))

colors = ["#E24B4A" if i == 0 else "#378ADD"
          for i in range(len(dist_trimmed))]

bars = ax.bar(range(len(dist_trimmed)), dist_trimmed.values,
              color=colors, alpha=0.85, edgecolor="white")

# Add value labels
for bar, val in zip(bars, dist_trimmed.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 20,
            f"{val:,}", ha="center", va="bottom", fontsize=9)

ax.set_xticks(range(len(dist_trimmed)))
ax.set_xticklabels([str(i) for i in dist_trimmed.index], fontsize=11)
ax.set_xlabel("Number of annotated isoforms in GENCODE v37", fontsize=13)
ax.set_ylabel("Number of always-100% dominant genes", fontsize=13)
ax.set_title("Isoform Count Distribution for Always-100% Dominant Genes\n"
             "(red = only 1 annotated isoform, blue = multiple isoforms)",
             fontsize=13)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
path2 = os.path.expanduser("~/finalproject/results/plots/plot_isoform_count_distribution.png")
plt.savefig(path2, dpi=150)
plt.close()
print(f"  Plot 2 saved to {path2}")

# ============================================
# SAVE RESULTS
# ============================================

always_100_with_counts.to_csv(
    os.path.expanduser("~/finalproject/results/always_100_genes.csv"), index=False)
print("\nFull list saved to ~/finalproject/results/always_100_genes.csv")
print("\nDone!")
