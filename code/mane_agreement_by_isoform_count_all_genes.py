import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ============================================
# CONFIGURATION
# ============================================
RESULTS_FILE = os.path.expanduser("~/finalproject/results/major_isoforms.csv")
GTF_FILE = "/opt/ccb/data/annotations/gencode.v37.annotation.gtf"

print("Loading results...")
results_df = pd.read_csv(RESULTS_FILE)
print(f"Loaded {len(results_df)} rows")

# ============================================
# STEP 1: COUNT ANNOTATED ISOFORMS FROM GTF
# ============================================

print("Parsing GTF for isoform counts...")

transcript_to_gene = {}

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

print(f"Done. {len(isoform_df)} protein-coding genes counted.")

# ============================================
# STEP 2: CALCULATE MANE AGREEMENT PER GENE
# ============================================

# Get one row per gene — majority MANE agreement across samples
gene_mane = (results_df
             .dropna(subset=["matches_mane"])
             .groupby("gene_id")["matches_mane"]
             .mean()
             .reset_index())
gene_mane.columns = ["gene_id", "mane_agreement"]

# Merge with isoform counts
merged = gene_mane.merge(isoform_df, on="gene_id", how="inner")

# Overall MANE agreement rate
overall_rate = merged["mane_agreement"].mean() * 100
print(f"Overall MANE agreement rate: {overall_rate:.1f}%")

# ============================================
# STEP 3: CALCULATE MANE AGREEMENT BY ISOFORM COUNT
# ============================================

# Group isoform counts — cap at 15+
merged["n_isoforms_grouped"] = merged["n_isoforms"].apply(
    lambda x: "15+" if x > 15 else str(x))

# Calculate mean MANE agreement per isoform count group
group_order = [str(i) for i in range(1, 16)] + ["15+"]
mane_by_isoform = (merged.groupby("n_isoforms_grouped")["mane_agreement"]
                   .agg(["mean", "count"])
                   .reindex(group_order)
                   .dropna())

# ============================================
# PLOT
# ============================================

fig, ax = plt.subplots(figsize=(14, 6))

bars = ax.bar(range(len(mane_by_isoform)),
              mane_by_isoform["mean"] * 100,
              color="#378ADD", alpha=0.85, edgecolor="white")

# Add value labels
for bar, val in zip(bars, mane_by_isoform["mean"] * 100):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{val:.1f}%", ha="center", va="bottom", fontsize=9)

# Add sample size labels below x axis
for i, (idx, row) in enumerate(mane_by_isoform.iterrows()):
    ax.text(i, -4, f"n={int(row['count'])}",
            ha="center", va="top", fontsize=8, color="gray")

# Overall MANE agreement line
ax.axhline(y=overall_rate, color="red", linestyle="--",
           linewidth=1.5, label=f"Overall MANE agreement ({overall_rate:.1f}%)")

ax.set_xticks(range(len(mane_by_isoform)))
ax.set_xticklabels(mane_by_isoform.index, fontsize=11)
ax.set_xlabel("Number of annotated isoforms in GENCODE v37", fontsize=13)
ax.set_ylabel("MANE_Select agreement rate (%)", fontsize=13)
ax.set_title("MANE_Select Agreement Rate for All Protein-Coding Genes\nby Annotated Isoform Count",
             fontsize=13)
ax.set_ylim(0, 105)
ax.legend(fontsize=11)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

output_path = os.path.expanduser(
    "~/finalproject/results/plots/plot_mane_agreement_by_isoform_count_all_genes.png")
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nPlot saved to {output_path}")
