"""
Remove proteins with C-terminal ER-retention motifs.
Motifs are supplied via snakemake.params.motifs.
Also checks dilysine motifs: KKXX and KXKX.
"""
import sys
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

MOTIFS = set(snakemake.params.motifs)


def is_er_retained(seq):
    if len(seq) < 4:
        return False
    tail = seq[-4:].upper()
    if tail in MOTIFS:
        return True
    # Dilysine KKXX
    if tail[0] == "K" and tail[1] == "K":
        return True
    # KXKX
    if tail[0] == "K" and tail[2] == "K":
        return True
    return False


df = pd.read_csv(snakemake.input.csv)

# Handle empty input
if df.empty or "Protein_ID" not in df.columns:
    df = pd.DataFrame(columns=["Protein_ID"])
    df.to_csv(snakemake.output.csv, index=False)
    open(snakemake.output.excluded, "w").close()
    open(snakemake.output.faa, "w").close()
    print("Input empty — wrote empty outputs", file=sys.stderr)
else:
    df["Protein_ID"] = df["Protein_ID"].astype(str)

    # Build seq dict
    seq_dict = {}
    records = list(SeqIO.parse(snakemake.input.faa, "fasta"))
    for rec in records:
        seq_dict[rec.id.split()[0]] = str(rec.seq)

    df["ER_Retained"] = df["Protein_ID"].map(lambda pid: is_er_retained(seq_dict.get(pid, "")))

    kept_df = df[~df["ER_Retained"]].drop(columns=["ER_Retained"]).copy()
    removed_df = df[df["ER_Retained"]].copy()

    kept_df.to_csv(snakemake.output.csv, index=False)
    with open(snakemake.output.excluded, "w") as fh:
        for pid in sorted(removed_df["Protein_ID"].unique()):
            fh.write(pid + "\n")

    keep_ids = set(kept_df["Protein_ID"])
    with open(snakemake.output.faa, "w") as fout:
        for rec in records:
            if rec.id.split()[0] in keep_ids:
                SeqIO.write(rec, fout, "fasta")

    print(f"Kept: {len(kept_df)} | ER removed: {len(removed_df)}", file=sys.stderr)
