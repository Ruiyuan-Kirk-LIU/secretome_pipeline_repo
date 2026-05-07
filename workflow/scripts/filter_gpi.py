"""
Remove GPI-anchored proteins predicted by NetGPI.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

# Parse NetGPI output
rows = []
with open(snakemake.input.gpi) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            rows.append({"Protein_ID": parts[0].split()[0],
                         "GPI_Prediction": parts[2].strip()})
gpi_df = pd.DataFrame(rows)

df = pd.read_csv(snakemake.input.csv)
df["Protein_ID"] = df["Protein_ID"].astype(str)
merged = df.merge(gpi_df, on="Protein_ID", how="left")
merged["GPI_Prediction"] = merged["GPI_Prediction"].fillna("NA")

keep_mask = merged["GPI_Prediction"].eq("Not GPI-Anchored")
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

print(f"Kept: {len(kept_df)} | GPI removed: {len(removed_df)}", file=sys.stderr)
