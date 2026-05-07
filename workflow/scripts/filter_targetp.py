"""
Remove proteins predicted by TargetP 2.0 to carry organelle-targeting signals.
Which prediction classes to remove is controlled by snakemake.params.remove_classes.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

REMOVE = set(snakemake.params.remove_classes)

# Parse TargetP output
rows = []
with open(snakemake.input.targetp) as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            rows.append({"Protein_ID": parts[0].split()[0],
                         "TargetP_Prediction": parts[1].strip()})
tp_df = pd.DataFrame(rows)

df = pd.read_csv(snakemake.input.csv)
df["Protein_ID"] = df["Protein_ID"].astype(str)
merged = df.merge(tp_df, on="Protein_ID", how="left")
merged["TargetP_Prediction"] = merged["TargetP_Prediction"].fillna("noTP")

keep_mask = ~merged["TargetP_Prediction"].isin(REMOVE)
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

print(f"Removed classes: {REMOVE} | Kept: {len(kept_df)} | Removed: {len(removed_df)}", file=sys.stderr)
