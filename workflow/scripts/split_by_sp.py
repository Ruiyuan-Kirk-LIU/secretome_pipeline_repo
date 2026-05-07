"""
Split nr99 proteome into SP+ and SP- based on SignalP 6.0 predictions.

SP+ = protein predicted as "SP" by SignalP 6.0.
SP- = all remaining proteins.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")


def read_signalp(path):
    """Return set of protein IDs predicted as SP by SignalP 6.0."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2 and parts[1].strip() == "SP":
                ids.add(parts[0].split()[0])
    return ids


# Collect SP-positive IDs
sp_ids = read_signalp(snakemake.input.signalp)

# Read all records
records = list(SeqIO.parse(snakemake.input.faa, "fasta"))

sp_pos, sp_neg = [], []
for rec in records:
    rid = rec.id.split()[0]
    if rid in sp_ids:
        sp_pos.append(rec)
    else:
        sp_neg.append(rec)

# Write FASTAs
SeqIO.write(sp_pos, snakemake.output.sp_pos_faa, "fasta")
SeqIO.write(sp_neg, snakemake.output.sp_neg_faa, "fasta")

# Write tracking CSVs
for label, recs, csv_path in [
    ("SP_pos", sp_pos, snakemake.output.sp_pos_csv),
    ("SP_neg", sp_neg, snakemake.output.sp_neg_csv),
]:
    df = pd.DataFrame({"Protein_ID": [r.id.split()[0] for r in recs]})
    df.to_csv(csv_path, index=False)

print(
    f"Total: {len(records)}  |  SP+: {len(sp_pos)}  |  SP-: {len(sp_neg)}",
    file=sys.stderr,
)
