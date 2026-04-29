import pandas as pd
import numpy as np
import os
import glob

# ============================================
# PART 1: PARSE GTF FILE
# ============================================

# Path to the GENCODE v37 GTF annotation file on the server
gtf_file = "/opt/ccb/data/annotations/gencode.v37.annotation.gtf"

# Set to store unique protein-coding gene IDs
# We use a set because we don't want duplicates
protein_coding_genes = set()

# Dictionary mapping each transcript ID to its parent gene ID
# e.g. {"ENST00000357654": "ENSG00000012048", ...}
transcript_to_gene = {}

# Dictionary mapping gene ID to its common name
# e.g. {"ENSG00000012048": "BRCA1", ...}
gene_to_name = {}

# Dictionary mapping each protein-coding gene to its MANE_Select transcript
# e.g. {"ENSG00000012048": "ENST00000357654", ...}
mane_select = {}

print("Parsing GTF file... this may take a minute")

# Open and read the GTF file line by line
# We read line by line instead of loading all at once because the file is ~1.5GB
with open(gtf_file) as f:
    for line in f:
        
        # Skip comment lines at the top of the GTF file
        # Comment lines start with "#" and contain metadata, not data
        if line.startswith("#"):
            continue
        
        # Split the line into its 9 columns using tab as separator
        fields = line.strip().split("\t")
        
        # Column 3 is the feature type (gene, transcript, exon, CDS etc.)
        feature_type = fields[2]
        
        # Column 9 contains all the attributes (gene_id, transcript_id, gene_type etc.)
        attributes = fields[8]
        
        # We only care about transcript lines
        # We skip gene, exon, CDS lines because transcript lines have all the info we need
        if feature_type != "transcript":
            continue
        
        # --- Extract gene_id from attributes ---
        # Attributes look like: gene_id "ENSG00000012048.23"; transcript_id "ENST00000357654.9"; ...
        gene_id = ""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1]
                # Remove version number after the dot
                # e.g. ENSG00000012048.23 → ENSG00000012048
                # We do this because Salmon quant.sf also strips version numbers
                gene_id = gene_id.split(".")[0]
        
        # --- Extract transcript_id from attributes ---
        transcript_id = ""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("transcript_id"):
                transcript_id = attr.split('"')[1]
                # Remove version number
                # e.g. ENST00000357654.9 → ENST00000357654
                transcript_id = transcript_id.split(".")[0]
        
        # --- Extract gene_name from attributes ---
        # gene_name is the human readable name like BRCA1, TP53 etc.
        gene_name = ""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1]
        
        # --- Extract gene_type from attributes ---
        # gene_type tells us if this is protein_coding, lncRNA, pseudogene etc.
        gene_type = ""
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_type"):
                gene_type = attr.split('"')[1]
        
        # --- Filter for protein-coding genes only ---
        # We only add genes and transcripts that are protein_coding
        # This is a key requirement of the project
        if gene_type == "protein_coding":
            # Add gene to our set of protein-coding genes
            protein_coding_genes.add(gene_id)
            # Map this transcript to its parent gene
            transcript_to_gene[transcript_id] = gene_id
            # Store the human readable gene name
            gene_to_name[gene_id] = gene_name
        
        # --- Check for MANE_Select tag ---
        # MANE_Select transcripts have the tag "MANE_Select" in their attributes
        # We only record MANE_Select for protein-coding genes
        if "MANE_Select" in attributes and gene_type == "protein_coding":
            mane_select[gene_id] = transcript_id

# Print summary to confirm everything worked
print(f"Done!")
print(f"Protein-coding genes found: {len(protein_coding_genes)}")
print(f"Transcripts mapped to genes: {len(transcript_to_gene)}")
print(f"MANE_Select transcripts found: {len(mane_select)}")


# ============================================
# PART 2: LOAD ALL 18 QUANT.SF FILES
# ============================================

# Path to the folder containing all sample quantification results
quants_dir = os.path.expanduser("~/finalproject/quants")

# Find all quant.sf files automatically using glob
# The * wildcard matches any sample name
quant_files = glob.glob(os.path.join(quants_dir, "*/quant.sf"))

# Sort the files so they are always in the same order
quant_files = sorted(quant_files)

print(f"\nFound {len(quant_files)} quant.sf files")

# Dictionary to store TPM data for each sample
# Key = sample name, Value = pandas Series of TPM values indexed by transcript ID
all_tpm = {}

print("Loading quant.sf files...")

for qfile in quant_files:
    
    # Extract the sample name from the file path
    # e.g. /home/sjameel1/finalproject/quants/ERR188023/quant.sf → ERR188023
    sample_name = os.path.basename(os.path.dirname(qfile))
    
    # Load the quant.sf file into a pandas dataframe
    # quant.sf is tab-separated with columns: Name, Length, EffectiveLength, TPM, NumReads
    df = pd.read_csv(qfile, sep="\t")
    
    # Clean up transcript IDs by removing version numbers
    # e.g. ENST00000357654.9 → ENST00000357654
    # This must match how we cleaned transcript IDs in the GTF
    df["Name"] = df["Name"].str.split(".").str[0]
    
    # Keep only the transcript name and TPM columns
    # Store as a Series indexed by transcript name
    all_tpm[sample_name] = df.set_index("Name")["TPM"]
    
    print(f"  Loaded {sample_name}: {len(df)} transcripts")

# Combine all samples into one big dataframe
# Rows = transcripts, Columns = samples
# This creates a matrix of TPM values
tpm_matrix = pd.DataFrame(all_tpm)

print(f"\nCombined TPM matrix shape: {tpm_matrix.shape}")
print(f"  → {tpm_matrix.shape[0]} transcripts x {tpm_matrix.shape[1]} samples")

# -----------------------------------------------
# Filter to protein-coding transcripts only
# -----------------------------------------------

# Keep only transcripts that are in our transcript_to_gene mapping
# (which only contains protein-coding transcripts from Part 1)
protein_coding_transcripts = set(transcript_to_gene.keys())

# Filter the matrix
tpm_matrix_pc = tpm_matrix[tpm_matrix.index.isin(protein_coding_transcripts)]

print(f"\nAfter filtering for protein-coding transcripts:")
print(f"  → {tpm_matrix_pc.shape[0]} transcripts x {tpm_matrix_pc.shape[1]} samples")

# Add a column mapping each transcript to its gene ID
# This uses our transcript_to_gene dictionary from Part 1
tpm_matrix_pc = tpm_matrix_pc.copy()
tpm_matrix_pc["gene_id"] = tpm_matrix_pc.index.map(transcript_to_gene)

print(f"\nPart 2 complete!")


# ============================================
# PART 3: FIND MAJOR ISOFORM PER GENE PER SAMPLE
# ============================================

print("\nCalculating major isoforms...")

# List to store results — one row per gene per sample
results = []

# Get list of all sample names (column names of our TPM matrix)
samples = [col for col in tpm_matrix_pc.columns if col != "gene_id"]

print(f"Total samples to process: {len(samples)}")

# Loop through each sample
for sample_num, sample in enumerate(samples, 1):
    
    # Progress tracker — shows which sample we are on
    print(f"\n[{sample_num}/{len(samples)}] Processing sample: {sample}...")
    
    # For this sample, get TPM values for all protein-coding transcripts
    # Create a dataframe with transcript ID, TPM value, and gene ID
    sample_df = pd.DataFrame({
        "transcript_id": tpm_matrix_pc.index,
        "TPM": tpm_matrix_pc[sample].values,
        "gene_id": tpm_matrix_pc["gene_id"].values
    })
    
    # Count how many genes we will process for this sample
    n_genes = sample_df["gene_id"].nunique()
    print(f"  → Analyzing {n_genes} protein-coding genes...")
    
    # Counter for expressed genes in this sample
    expressed_genes = 0
    
    # Group transcripts by gene
    # For each gene, we want to know:
    # 1. Total TPM (sum of all isoforms)
    # 2. Which transcript has the highest TPM (major isoform)
    # 3. What is that highest TPM value
    for gene_id, gene_transcripts in sample_df.groupby("gene_id"):
        
        # Calculate total expression at this gene locus
        total_tpm = gene_transcripts["TPM"].sum()
        
        # Skip genes with zero total expression
        # These genes are not expressed in this sample at all
        if total_tpm == 0:
            continue
        
        # Count expressed genes
        expressed_genes += 1
        
        # Find the major isoform — transcript with highest TPM
        major_idx = gene_transcripts["TPM"].idxmax()
        major_transcript = gene_transcripts.loc[major_idx, "transcript_id"]
        major_tpm = gene_transcripts.loc[major_idx, "TPM"]
        
        # Calculate percentage contribution of major isoform
        major_pct = (major_tpm / total_tpm) * 100
        
        # Is this isoform dominant? (contributes more than 50%)
        is_dominant = major_pct > 50
        
        # Get the human readable gene name
        gene_name = gene_to_name.get(gene_id, "unknown")
        
        # Get the MANE_Select transcript for this gene (if it exists)
        mane_transcript = mane_select.get(gene_id, None)
        
        # Does the major isoform match the MANE_Select transcript?
        matches_mane = (major_transcript == mane_transcript) if mane_transcript else None
        
        # Store this result
        results.append({
            "sample": sample,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "total_tpm": total_tpm,
            "major_transcript": major_transcript,
            "major_tpm": major_tpm,
            "major_pct": major_pct,
            "is_dominant": is_dominant,
            "mane_transcript": mane_transcript,
            "matches_mane": matches_mane
        })
    
    # Summary for this sample
    print(f"  → Expressed genes: {expressed_genes} out of {n_genes}")
    print(f"  → Sample {sample} done! ✓")

# Convert results list to a pandas dataframe
results_df = pd.DataFrame(results)

print(f"\n{'='*50}")
print(f"ALL SAMPLES PROCESSED!")
print(f"{'='*50}")
print(f"Total gene-sample combinations analyzed: {len(results_df)}")
print(f"\nFirst few rows of results:")
print(results_df.head(10).to_string())

# Save results to a CSV file for later use
results_df.to_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"), index=False)
print(f"\nResults saved to ~/finalproject/results/major_isoforms.csv")
print(f"\nPart 3 complete!")


# ============================================
# PART 4: ANSWER THE 4 RESEARCH QUESTIONS
# ============================================

print("\n" + "="*50)
print("PART 4: RESEARCH QUESTIONS")
print("="*50)

# -----------------------------------------------
# RESEARCH QUESTION 1:
# Do most protein-coding genes have a single dominant isoform?
# -----------------------------------------------

print("\n--- Research Question 1 ---")
print("Do most protein-coding genes have a dominant isoform (>50%)?")

# For each sample, count how many genes have a dominant isoform
q1_per_sample = results_df.groupby("sample")["is_dominant"].agg(
    total_genes="count",
    dominant_genes="sum"
)

# Calculate percentage of dominant genes per sample
q1_per_sample["pct_dominant"] = (q1_per_sample["dominant_genes"] / q1_per_sample["total_genes"]) * 100

print("\nPer sample summary:")
print(q1_per_sample.to_string())

# Overall average across all samples
avg_pct_dominant = q1_per_sample["pct_dominant"].mean()
print(f"\nOn average across all samples:")
print(f"  {avg_pct_dominant:.1f}% of protein-coding genes have a dominant isoform (>50%)")

# Distribution of major isoform percentages
print(f"\nDistribution of major isoform % across all genes and samples:")
print(results_df["major_pct"].describe())

# Save Q1 results
q1_per_sample.to_csv(os.path.expanduser("~/finalproject/results/q1_dominance_per_sample.csv"))
print("\nQ1 results saved!")

# -----------------------------------------------
# RESEARCH QUESTION 2:
# How stable is the major isoform identity across individuals?
# -----------------------------------------------

print("\n--- Research Question 2 ---")
print("How stable is the major isoform identity across samples?")

# For each gene, find how many unique transcripts are called major isoform
# across all 18 samples
q2 = results_df.groupby("gene_id").agg(
    gene_name=("gene_name", "first"),
    n_samples=("sample", "count"),
    n_unique_major=("major_transcript", "nunique")
).reset_index()

# A gene is consistent if it always uses the same major isoform
# i.e. n_unique_major == 1
q2["is_consistent"] = q2["n_unique_major"] == 1

# Summary
n_consistent = q2["is_consistent"].sum()
n_total = len(q2)
pct_consistent = (n_consistent / n_total) * 100

print(f"\nOut of {n_total} expressed protein-coding genes:")
print(f"  {n_consistent} ({pct_consistent:.1f}%) always use the same major isoform across all samples")
print(f"  {n_total - n_consistent} ({100-pct_consistent:.1f}%) switch their major isoform at least once")

# Distribution of number of unique major isoforms per gene
print(f"\nDistribution of unique major isoforms per gene:")
print(q2["n_unique_major"].value_counts().sort_index().to_string())

# Save Q2 results
q2.to_csv(os.path.expanduser("~/finalproject/results/q2_isoform_consistency.csv"), index=False)
print("\nQ2 results saved!")

# -----------------------------------------------
# RESEARCH QUESTION 3:
# What proportion of dominant isoforms correspond to MANE_Select?
# -----------------------------------------------

print("\n--- Research Question 3 ---")
print("What proportion of dominant isoforms match MANE_Select?")

# Filter to only genes that:
# 1. Have a dominant isoform (>50%)
# 2. Have a MANE_Select transcript (not None)
dominant_with_mane = results_df[
    (results_df["is_dominant"] == True) &
    (results_df["matches_mane"].notna())
]

# Count matches
n_matches = dominant_with_mane["matches_mane"].sum()
n_total_dominant = len(dominant_with_mane)
pct_matches = (n_matches / n_total_dominant) * 100

print(f"\nAmong dominant isoforms with a MANE_Select transcript:")
print(f"  Total dominant gene-sample pairs: {n_total_dominant}")
print(f"  Match MANE_Select: {n_matches} ({pct_matches:.1f}%)")
print(f"  Differ from MANE_Select: {n_total_dominant - n_matches} ({100-pct_matches:.1f}%)")

# Also look at non-dominant genes
non_dominant_with_mane = results_df[
    (results_df["is_dominant"] == False) &
    (results_df["matches_mane"].notna())
]
n_matches_nd = non_dominant_with_mane["matches_mane"].sum()
n_total_nd = len(non_dominant_with_mane)
pct_matches_nd = (n_matches_nd / n_total_nd) * 100

print(f"\nAmong NON-dominant isoforms with a MANE_Select transcript:")
print(f"  Total non-dominant gene-sample pairs: {n_total_nd}")
print(f"  Match MANE_Select: {n_matches_nd} ({pct_matches_nd:.1f}%)")

# Save Q3 results
q3_summary = pd.DataFrame({
    "category": ["dominant", "non_dominant"],
    "total": [n_total_dominant, n_total_nd],
    "matches_mane": [n_matches, n_matches_nd],
    "pct_matches": [pct_matches, pct_matches_nd]
})
q3_summary.to_csv(os.path.expanduser("~/finalproject/results/q3_mane_comparison.csv"), index=False)
print("\nQ3 results saved!")

# -----------------------------------------------
# RESEARCH QUESTION 4:
# Are there genes where the major isoform frequently switches?
# -----------------------------------------------

print("\n--- Research Question 4 ---")
print("Which genes most frequently switch their major isoform?")

# Use Q2 results — genes with most unique major isoforms are the most switching
# Sort by number of unique major isoforms (descending)
switching_genes = q2[q2["n_unique_major"] > 1].sort_values(
    "n_unique_major", ascending=False
)

print(f"\nTotal switching genes: {len(switching_genes)}")
print(f"\nTop 20 most frequently switching genes:")
print(switching_genes[["gene_name", "n_samples", "n_unique_major"]].head(20).to_string())

# Save Q4 results
switching_genes.to_csv(os.path.expanduser("~/finalproject/results/q4_switching_genes.csv"), index=False)
print("\nQ4 results saved!")

print("\n" + "="*50)
print("ALL RESEARCH QUESTIONS ANSWERED!")
print("Results saved to ~/finalproject/results/")
print("="*50)


# ============================================
# PART 5: PLOTTING
# ============================================

import matplotlib
# Use non-interactive backend since server has no display
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Create plots directory
plots_dir = os.path.expanduser("~/finalproject/results/plots")
os.makedirs(plots_dir, exist_ok=True)

# Set a clean style for all plots
sns.set_style("whitegrid")
sns.set_palette("husl")

print("\n" + "="*50)
print("PART 5: PLOTTING")
print("="*50)

# -----------------------------------------------
# PLOT 1: Distribution of major isoform %
# -----------------------------------------------

print("\nPlotting Plot 1: Distribution of major isoform %...")

fig, ax = plt.subplots(figsize=(10, 6))

# Plot histogram of major isoform % across all genes and samples
ax.hist(results_df["major_pct"], bins=50, color="steelblue", edgecolor="white", alpha=0.8)

# Add vertical line at 50% threshold
ax.axvline(x=50, color="red", linestyle="--", linewidth=2, label="50% dominance threshold")

# Labels and title
ax.set_xlabel("Major Isoform % of Total Gene Expression", fontsize=13)
ax.set_ylabel("Number of Gene-Sample Combinations", fontsize=13)
ax.set_title("Distribution of Major Isoform Expression Percentage\nAcross All Protein-Coding Genes and Samples", fontsize=14)
ax.legend(fontsize=12)

# Add text showing % above and below threshold
pct_above = (results_df["major_pct"] > 50).mean() * 100
ax.text(0.05, 0.95, f"{pct_above:.1f}% of gene-sample pairs\nhave dominant isoform (>50%)",
        transform=ax.transAxes, fontsize=11, verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot1_major_isoform_distribution.png"), dpi=150)
plt.close()
print("  Plot 1 saved!")

# -----------------------------------------------
# PLOT 2: % of genes with dominant isoform per sample
# -----------------------------------------------

print("Plotting Plot 2: Dominance per sample...")

fig, ax = plt.subplots(figsize=(12, 6))

# Plot bar chart of % dominant genes per sample
bars = ax.bar(q1_per_sample.index, q1_per_sample["pct_dominant"],
              color="steelblue", edgecolor="white", alpha=0.8)

# Add value labels on top of each bar
for bar, val in zip(bars, q1_per_sample["pct_dominant"]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{val:.1f}%", ha="center", va="bottom", fontsize=9)

# Add horizontal line showing average
avg = q1_per_sample["pct_dominant"].mean()
ax.axhline(y=avg, color="red", linestyle="--", linewidth=2,
           label=f"Average: {avg:.1f}%")

# Labels and title
ax.set_xlabel("Sample", fontsize=13)
ax.set_ylabel("% of Genes with Dominant Isoform (>50%)", fontsize=13)
ax.set_title("Percentage of Protein-Coding Genes with a Dominant Isoform\nPer Sample", fontsize=14)
ax.set_xticklabels(q1_per_sample.index, rotation=45, ha="right", fontsize=9)
ax.legend(fontsize=12)
ax.set_ylim(0, 100)

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot2_dominance_per_sample.png"), dpi=150)
plt.close()
print("  Plot 2 saved!")

# -----------------------------------------------
# PLOT 3: Consistent vs switching genes
# -----------------------------------------------

print("Plotting Plot 3: Consistent vs switching genes...")

fig, ax = plt.subplots(figsize=(8, 8))

# Count consistent vs switching genes
n_consistent = q2["is_consistent"].sum()
n_switching = (~q2["is_consistent"]).sum()

# Pie chart
wedges, texts, autotexts = ax.pie(
    [n_consistent, n_switching],
    labels=[f"Consistent\n({n_consistent:,} genes)", f"Switching\n({n_switching:,} genes)"],
    autopct="%1.1f%%",
    colors=["steelblue", "salmon"],
    startangle=90,
    textprops={"fontsize": 12}
)

ax.set_title("Stability of Major Isoform Identity\nAcross All 18 Samples", fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot3_consistent_vs_switching.png"), dpi=150)
plt.close()
print("  Plot 3 saved!")

# -----------------------------------------------
# PLOT 4: Distribution of unique major isoforms per gene
# -----------------------------------------------

print("Plotting Plot 4: Distribution of unique major isoforms per gene...")

fig, ax = plt.subplots(figsize=(10, 6))

# Count genes with each number of unique major isoforms
unique_counts = q2["n_unique_major"].value_counts().sort_index()

bars = ax.bar(unique_counts.index, unique_counts.values,
              color="steelblue", edgecolor="white", alpha=0.8)

# Add value labels on top of each bar
for bar, val in zip(bars, unique_counts.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
            f"{val:,}", ha="center", va="bottom", fontsize=9)

# Labels and title
ax.set_xlabel("Number of Unique Major Isoforms Across 18 Samples", fontsize=13)
ax.set_ylabel("Number of Genes", fontsize=13)
ax.set_title("Distribution of Major Isoform Switching\nAcross Protein-Coding Genes", fontsize=14)
ax.set_xticks(unique_counts.index)

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot4_unique_isoforms_distribution.png"), dpi=150)
plt.close()
print("  Plot 4 saved!")

# -----------------------------------------------
# PLOT 5: MANE_Select comparison
# -----------------------------------------------

print("Plotting Plot 5: MANE_Select comparison...")

fig, ax = plt.subplots(figsize=(8, 6))

# Compare dominant vs non-dominant genes for MANE matching
categories = ["Dominant\nIsoforms\n(>50%)", "Non-Dominant\nIsoforms\n(<50%)"]
match_pcts = [pct_matches, pct_matches_nd]
no_match_pcts = [100 - pct_matches, 100 - pct_matches_nd]

x = range(len(categories))
width = 0.35

bars1 = ax.bar([i - width/2 for i in x], match_pcts, width,
               label="Matches MANE_Select", color="steelblue", alpha=0.8)
bars2 = ax.bar([i + width/2 for i in x], no_match_pcts, width,
               label="Differs from MANE_Select", color="salmon", alpha=0.8)

# Add value labels
for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=11)
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
            f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=11)

ax.set_ylabel("Percentage (%)", fontsize=13)
ax.set_title("Comparison of Major Isoforms to MANE_Select Transcripts", fontsize=14)
ax.set_xticks(list(x))
ax.set_xticklabels(categories, fontsize=12)
ax.legend(fontsize=11)
ax.set_ylim(0, 110)

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot5_mane_comparison.png"), dpi=150)
plt.close()
print("  Plot 5 saved!")

# -----------------------------------------------
# PLOT 6: Top 20 most frequently switching genes
# -----------------------------------------------

print("Plotting Plot 6: Top 20 switching genes...")

fig, ax = plt.subplots(figsize=(12, 7))

# Get top 20 switching genes
top20 = switching_genes.head(20)

bars = ax.barh(top20["gene_name"], top20["n_unique_major"],
               color="salmon", edgecolor="white", alpha=0.8)

# Add value labels
for bar, val in zip(bars, top20["n_unique_major"]):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
            f"{val}", ha="left", va="center", fontsize=10)

ax.set_xlabel("Number of Unique Major Isoforms Across 18 Samples", fontsize=13)
ax.set_ylabel("Gene", fontsize=13)
ax.set_title("Top 20 Genes with Most Frequent Major Isoform Switching", fontsize=14)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(os.path.join(plots_dir, "plot6_top_switching_genes.png"), dpi=150)
plt.close()
print("  Plot 6 saved!")

print("\n" + "="*50)
print("ALL PLOTS SAVED!")
print(f"Location: {plots_dir}")
print("="*50)
