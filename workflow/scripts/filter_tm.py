"""
Remove proteins with predicted transmembrane helices (DeepTMHMM).
Proteins whose sole TM helix overlaps the N-terminal signal peptide region
(first 70 aa) are exempt and retained.
"""
import sys, re
import pandas as pd
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

SP_REGION = 70  # max length of signal peptide region (aa)

# --- Parse DeepTMHMM GFF3 ---
tm_count = {}
tm_coords = {}  # {pid: [(start, end), ...]}
pat = re.compile(r"#\s+(\S+)\s+Number of predicted TMRs:\s+(\d+)")

with open(snakemake.input.gff3) as fh:
    for line in fh:
        m = pat.match(line.strip())
        if m:
            pid, n = m.group(1), int(m.group(2))
            tm_count[pid] = n
        elif not line.startswith("#"):
            cols = line.strip().split("\t")
            if len(cols) >= 5 and cols[2] == "TMhelix":
                pid = cols[0]
                tm_coords.setdefault(pid, []).append((int(cols[3]), int(cols[4])))

# --- Determine which proteins to keep ---
df = pd.read_csv(snakemake.input.csv)
df["Protein_ID"] = df["Protein_ID"].astype(str)
df["Predicted_TMRs"] = df["Protein_ID"].map(tm_count).fillna(0).astype(int)

def is_sp_only_tm(pid):
    """Return True if protein has exactly 1 TM helix within the SP region."""
    coords = tm_coords.get(pid, [])
    if len(coords) != 1:
        return False
    start, end = coords[0]
    return start <= SP_REGION and end <= SP_REGION

keep_mask = (df["Predicted_TMRs"] == 0) | df["Protein_ID"].apply(is_sp_only_tm)
kept_df = df[keep_mask].copy()
removed_df = df[~keep_mask].copy()

kept_df.to_csv(snakemake.output.csv, index=False)

# Write excluded IDs
with open(snakemake.output.excluded, "w") as fh:
    for pid in sorted(removed_df["Protein_ID"].unique()):
        fh.write(pid + "\n")

# Write FASTA subset
keep_ids = set(kept_df["Protein_ID"])
n = 0
with open(snakemake.output.faa, "w") as fout:
    for rec in SeqIO.parse(snakemake.input.faa, "fasta"):
        if rec.id.split()[0] in keep_ids:
            SeqIO.write(rec, fout, "fasta")
            n += 1

print(f"Input: {len(df)} | Kept: {len(kept_df)} | TM removed: {len(removed_df)}", file=sys.stderr)
