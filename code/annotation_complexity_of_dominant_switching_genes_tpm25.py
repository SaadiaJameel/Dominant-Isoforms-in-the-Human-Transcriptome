import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# ============================================
# CONFIGURATION
# ============================================
RESULTS_FILE = os.path.expanduser("~/finalproject/results/major_isoforms_tpm25.csv")
GTF_FILE = "/opt/ccb/data/annotations/gencode.v37.annotation.gtf"
MIN_TPM = 25

# ============================================
# STEP 1: COUNT ANNOTATED ISOFORMS FROM GTF
# ============================================

print("Counting annotated isoforms from GTF...")

transcript_to_gene = {}
gene_to_name = {}

with open(GTF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "transcript":
            continue
        attributes = fields[8]

        gene_id = ""
        transcript_id = ""
        gene_name = ""
        gene_type = ""

        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1].split(".")[0]
            if attr.startswith("transcript_id"):
                transcript_id = attr.split('"')[1].split(".")[0]
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1]
            if attr.startswith("gene_type"):
                gene_type = attr.split('"')[1]

        if gene_type == "protein_coding":
            transcript_to_gene[transcript_id] = gene_id
            gene_to_name[gene_id] = gene_name

print(f"Done. {len(transcript_to_gene)} protein-coding transcripts counted.")

isoform_counts = (pd.Series(transcript_to_gene)
                  .reset_index()
                  .rename(columns={"index": "transcript_id", 0: "gene_id"})
                  .groupby("gene_id")
                  .size()
                  .reset_index(name="n_annotated_isoforms"))

isoform_counts["gene_name"] = isoform_counts["gene_id"].map(gene_to_name)

# ============================================
# STEP 2: LOAD RESULTS AND FILTER TO DOMINANT
# ============================================

print("Loading results...")
results_df = pd.read_csv(RESULTS_FILE)
results_df = results_df[results_df["is_dominant"] == True]
print(f"After filtering to dominant isoforms: {len(results_df)} gene-sample pairs")

q2 = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

# ============================================
# STEP 3: TOP 20 SWITCHING, BOTTOM 20 STABLE
# ============================================

top20 = q2[q2["n_unique_major"] > 1].sort_values(
    "n_unique_major", ascending=False).head(20).copy()

bottom20 = q2[q2["n_unique_major"] == 1].sample(
    20, random_state=42).copy()

# Merge with annotated isoform counts
top20 = top20.merge(isoform_counts[["gene_id", "n_annotated_isoforms"]], on="gene_id", how="left")
bottom20 = bottom20.merge(isoform_counts[["gene_id", "n_annotated_isoforms"]], on="gene_id", how="left")

top20["group"] = "Top 20 (Most Switching)"
bottom20["group"] = "Bottom 20 (Stable)"

combined = pd.concat([top20, bottom20], ignore_index=True)

print(f"\nTotal genes in plot: {len(combined)}")
print("\nTop 20:")
print(top20[["gene_name", "n_unique_major", "n_annotated_isoforms"]].to_string())
print("\nBottom 20 (stable):")
print(bottom20[["gene_name", "n_unique_major", "n_annotated_isoforms"]].to_string())

# ============================================
# PLOT
# ============================================

colors = ["#E24B4A" if g == "Top 20 (Most Switching)" else "#378ADD"
          for g in combined["group"]]

fig, ax = plt.subplots(figsize=(12, 14))

bars = ax.barh(range(len(combined)), combined["n_annotated_isoforms"],
               color=colors, edgecolor="white", alpha=0.85)

for i, (bar, val) in enumerate(zip(bars, combined["n_annotated_isoforms"])):
    ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
            f"{int(val)}", ha="left", va="center", fontsize=9)

ax.set_yticks(range(len(combined)))
ax.set_yticklabels(combined["gene_name"], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel("Number of Annotated Isoforms in GENCODE v37", fontsize=13)
ax.set_ylabel("Gene", fontsize=13)
ax.set_title(f"Annotated Isoform Count:\nTop 20 Switching vs 20 Stable Genes (TPM >= {MIN_TPM}, Dominant Only)",
             fontsize=14)

legend_elements = [Patch(facecolor="#E24B4A", label="Top 20 (Most Switching)"),
                   Patch(facecolor="#378ADD", label="Stable Genes")]
ax.legend(handles=legend_elements, fontsize=11, loc="lower right")

plt.tight_layout()

output_path = os.path.expanduser(
    f"~/finalproject/results/plots/plot_top_bottom20_annotated_isoforms_dominant_tpm{MIN_TPM}.png")
plt.savefig(output_path, dpi=150)
plt.close()
print(f"\nPlot saved to {output_path}")
