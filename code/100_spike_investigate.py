import pandas as pd
import numpy as np
import os

# Load the already-computed results — no need to rerun Salmon or parse GTF!
print("Loading results...")
results_df = pd.read_csv(os.path.expanduser("~/finalproject/results/major_isoforms.csv"))
print(f"Loaded {len(results_df)} rows")

# ============================================
# INVESTIGATING THE 100% SPIKE
# ============================================

print("\n--- Investigating the 100% dominance spike ---")

genes_100 = results_df[results_df["major_pct"] == 100]
genes_partial = results_df[(results_df["major_pct"] > 50) & (results_df["major_pct"] < 100)]
genes_non_dominant = results_df[results_df["major_pct"] <= 50]

print(f"100% dominant gene-sample pairs:       {len(genes_100)}")
print(f"Partial dominant (50-100%) pairs:      {len(genes_partial)}")
print(f"Non-dominant (<50%) pairs:             {len(genes_non_dominant)}")

# MANE agreement rates
genes_100_mane = genes_100.dropna(subset=["matches_mane"])
genes_partial_mane = genes_partial.dropna(subset=["matches_mane"])

mane_100 = genes_100_mane["matches_mane"].mean() * 100
mane_partial = genes_partial_mane["matches_mane"].mean() * 100

print(f"\nMANE_Select agreement rate:")
print(f"  100% dominant genes:      {mane_100:.1f}%")
print(f"  Partially dominant genes: {mane_partial:.1f}%")

if mane_100 > mane_partial:
    print("  → 100% genes agree with MANE more → likely REAL BIOLOGY")
else:
    print("  → 100% genes agree with MANE less → possibly ARTIFACT")

# Consistency across samples
print("\n--- Consistency of 100% across samples ---")
pct_100_consistency = results_df.groupby("gene_id").apply(
    lambda x: (x["major_pct"] == 100).mean() * 100
).reset_index()
pct_100_consistency.columns = ["gene_id", "pct_samples_at_100"]

always_100 = (pct_100_consistency["pct_samples_at_100"] == 100).sum()
sometimes_100 = ((pct_100_consistency["pct_samples_at_100"] > 0) &
                 (pct_100_consistency["pct_samples_at_100"] < 100)).sum()
never_100 = (pct_100_consistency["pct_samples_at_100"] == 0).sum()

print(f"  Always 100% across all 18 samples:  {always_100} genes")
print(f"  100% in only some samples:          {sometimes_100} genes")
print(f"  Never 100% in any sample:           {never_100} genes")

# Save
pct_100_consistency.to_csv(
    os.path.expanduser("~/finalproject/results/bonus_100pct_analysis.csv"), index=False)
print("\nResults saved to ~/finalproject/results/bonus_100pct_analysis.csv")
