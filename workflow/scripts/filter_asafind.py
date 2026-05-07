"""
Remove proteins predicted by ASAFind to be plastid-targeted (complex plastid organisms).
Retains proteins predicted as 'Not plastid' by ASAFind 2.0.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

# Parse ASAFind output (tab-separated)
asafind = pd.read_csv(snakemake.input.asafind, sep="\t")
asafind = asafind.rename(columns={"Identifier": "Protein_ID"})

# Find the prediction column (contains "Prediction" to avoid matching
# "ASAFind cleavage position" or "ASAFind 20aa score")
pred_col = next((c for c in asafind.columns if "Prediction" in c), None)
if pred_col is None:
    raise ValueError(f"Cannot find ASAFind prediction column in: {list(asafind.columns)}")
asafind = asafind[["Protein_ID", pred_col]].rename(columns={pred_col: "ASAFind_Prediction"})

print(f"ASAFind prediction column: {pred_col}", file=sys.stderr)
print(f"Unique predictions: {asafind['ASAFind_Prediction'].unique()}", file=sys.stderr)

df = pd.read_csv(snakemake.input.csv)
df["Protein_ID"] = df["Protein_ID"].astype(str)
merged = df.merge(asafind, on="Protein_ID", how="left")
merged["ASAFind_Prediction"] = merged["ASAFind_Prediction"].fillna("NA")

# Keep proteins predicted as "Not plastid" (covers both ASAFind 1.x and 2.0 labels)
keep_mask = merged["ASAFind_Prediction"].str.contains("Not plastid", case=False, na=False)
kept_df = merged[keep_mask].copy()
removed_df = merged[~keep_mask].copy()

kept_df.to_csv(snakemake.output.csv, index=False)
with open(snakemake.output.excluded, "w") as fh:
    for pid in sorted(removed_df["Protein_ID"].unique()):
        fh.write(pid + "\n")

keep_ids = set(kept_df["Protein_ID"])
with open(snakemake.output.faa, "w") as fout:
    for rec in SeqIO.parse(snakemake.input.faa, "fasta"):
        if rec.id.split()[0] in keep_ids:
            SeqIO.write(rec, fout, "fasta")

print(f"Kept: {len(kept_df)} | Plastid removed: {len(removed_df)}", file=sys.stderr)
