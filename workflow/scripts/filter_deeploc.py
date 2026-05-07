"""
Retain only proteins predicted as extracellular by DeepLoc 2.0.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

deeploc = pd.read_csv(snakemake.input.deeploc)
extra_ids = set(
    deeploc.loc[
        deeploc["Localizations"].astype(str).str.contains("Extracellular", case=False, na=False),
        "Protein_ID",
    ].astype(str)
)

df = pd.read_csv(snakemake.input.csv)
df["Protein_ID"] = df["Protein_ID"].astype(str)
kept_df = df[df["Protein_ID"].isin(extra_ids)].copy()
kept_df.to_csv(snakemake.output.csv, index=False)

keep_ids = set(kept_df["Protein_ID"])
n = 0
with open(snakemake.output.faa, "w") as fout:
    for rec in SeqIO.parse(snakemake.input.faa, "fasta"):
        if rec.id.split()[0] in keep_ids:
            SeqIO.write(rec, fout, "fasta")
            n += 1

print(f"Extracellular: {len(kept_df)} / {len(df)} | FASTA: {n}", file=sys.stderr)
