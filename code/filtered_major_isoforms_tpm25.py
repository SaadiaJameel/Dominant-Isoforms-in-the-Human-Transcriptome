import pandas as pd
import numpy as np
import os
import glob

# ============================================
# CONFIGURATION
# ============================================

MIN_TPM = 25  # minimum total gene TPM — change to 5 or 10 to test other thresholds

gtf_file = "/opt/ccb/data/annotations/gencode.v37.annotation.gtf"
quants_dir = os.path.expanduser("~/finalproject/quants")
output_file = os.path.expanduser(f"~/finalproject/results/major_isoforms_tpm{MIN_TPM}.csv")

# ============================================
# PART 1: PARSE GTF FILE
# ============================================

print("Parsing GTF file... this may take a minute")

protein_coding_genes = set()
transcript_to_gene = {}
gene_to_name = {}
mane_select = {}

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
            protein_coding_genes.add(gene_id)
            transcript_to_gene[transcript_id] = gene_id
            gene_to_name[gene_id] = gene_name

        if "MANE_Select" in attributes and gene_type == "protein_coding":
            mane_select[gene_id] = transcript_id

print(f"Protein-coding genes:      {len(protein_coding_genes)}")
print(f"Transcripts mapped:        {len(transcript_to_gene)}")
print(f"MANE_Select transcripts:   {len(mane_select)}")

# ============================================
# PART 2: LOAD ALL 18 QUANT.SF FILES
# ============================================

print("\nLoading quant.sf files...")

quant_files = sorted(glob.glob(os.path.join(quants_dir, "*/quant.sf")))
print(f"Found {len(quant_files)} quant.sf files")

all_tpm = {}
protein_coding_transcripts = set(transcript_to_gene.keys())

for qfile in quant_files:
    sample_name = os.path.basename(os.path.dirname(qfile))
    df = pd.read_csv(qfile, sep="\t")

    # Clean transcript IDs — remove version numbers
    df["Name"] = df["Name"].str.split(".").str[0]

    # Filter to protein-coding transcripts only
    df = df[df["Name"].isin(protein_coding_transcripts)]

    all_tpm[sample_name] = df.set_index("Name")["TPM"]
    print(f"  Loaded {sample_name}: {len(df)} protein-coding transcripts")

# Combine into one matrix
tpm_matrix = pd.DataFrame(all_tpm)
tpm_matrix["gene_id"] = tpm_matrix.index.map(transcript_to_gene)

print(f"\nTPM matrix shape: {tpm_matrix.shape}")

# ============================================
# PART 3: CALCULATE MAJOR ISOFORM
# WITH TPM FILTER APPLIED BEFORE CALCULATION
# ============================================

print(f"\nCalculating major isoforms (TPM >= {MIN_TPM} filter applied)...")

samples = [col for col in tpm_matrix.columns if col != "gene_id"]
results = []

for sample_num, sample in enumerate(samples, 1):
    print(f"  [{sample_num}/{len(samples)}] Processing {sample}...")

    sample_df = pd.DataFrame({
        "transcript_id": tpm_matrix.index,
        "TPM": tpm_matrix[sample].values,
        "gene_id": tpm_matrix["gene_id"].values
    })

    expressed_genes = 0
    filtered_genes = 0

    for gene_id, gene_transcripts in sample_df.groupby("gene_id"):

        # Sum all isoform TPMs for this gene
        total_tpm = gene_transcripts["TPM"].sum()

        # ── KEY FILTER ──────────────────────────────
        # Skip genes with total expression below threshold
        # This is applied BEFORE major isoform calculation
        # so lowly expressed genes never enter the analysis
        if total_tpm < MIN_TPM:
            filtered_genes += 1
            continue
        # ────────────────────────────────────────────

        expressed_genes += 1

        # Find major isoform
        major_idx = gene_transcripts["TPM"].idxmax()
        major_transcript = gene_transcripts.loc[major_idx, "transcript_id"]
        major_tpm = gene_transcripts.loc[major_idx, "TPM"]

        # Calculate percentage
        major_pct = (major_tpm / total_tpm) * 100

        # Is it dominant?
        is_dominant = major_pct > 50

        # Gene name and MANE info
        gene_name = gene_to_name.get(gene_id, "unknown")
        mane_transcript = mane_select.get(gene_id, None)
        matches_mane = (major_transcript == mane_transcript) if mane_transcript else None

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

    print(f"    Expressed (>= {MIN_TPM} TPM): {expressed_genes} genes")
    print(f"    Filtered out (< {MIN_TPM} TPM): {filtered_genes} genes")

# Convert to dataframe
results_df = pd.DataFrame(results)

print(f"\nTotal gene-sample pairs: {len(results_df)}")
print(f"Unique genes:            {results_df['gene_id'].nunique()}")

# ============================================
# SAVE
# ============================================

results_df.to_csv(output_file, index=False)
print(f"\nSaved to {output_file}")
print("Done!")
