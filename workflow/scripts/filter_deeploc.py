"""
Retain only proteins predicted as extracellular.
Supports both DeepLoc 2.1 (eukaryotes) and DeepLocPro 1.0 (prokaryotes).
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

deeploc = pd.read_csv(snakemake.input.deeploc)

# Find the protein ID column (try common names)
id_col = None
for candidate in ["Protein_ID", "protein_id", "ID", "id", "Protein ID"]:
    if candidate in deeploc.columns:
        id_col = candidate
        break
if id_col is None:
    id_col = deeploc.columns[0]
    print(f"WARNING: No standard ID column found, using first column: {id_col}", file=sys.stderr)

# Find the localisation column (try common names from DeepLoc 2.1 and DeepLocPro 1.0)
loc_col = None
for candidate in ["Localizations", "Localization", "Location", "Prediction",
                   "localizations", "localization", "location", "prediction"]:
    if candidate in deeploc.columns:
        loc_col = candidate
        break
if loc_col is None:
    # Check for any column containing "local" or "predict" in its name
    for col in deeploc.columns:
        if "local" in col.lower() or "predict" in col.lower():
            loc_col = col
            break
if loc_col is None:
    raise ValueError(f"Cannot find localisation column in: {list(deeploc.columns)}")

print(f"Using ID column: {id_col}, Location column: {loc_col}", file=sys.stderr)

# Extract extracellular protein IDs
extra_ids = set(
    deeploc.loc[
        deeploc[loc_col].astype(str).str.contains("Extracellular", case=False, na=False),
        id_col,
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
