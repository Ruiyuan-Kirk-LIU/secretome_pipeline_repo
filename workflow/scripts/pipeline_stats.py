"""
Collect pipeline summary statistics for one sample.
"""
import sys
import os
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w")

sample = snakemake.wildcards.sample
organism_type = snakemake.config["samples"][sample]["organism_type"]


def count_fasta(path):
    if not os.path.exists(path):
        return "NA"
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


base = f"results/{sample}"
rows = []

rows.append(("nr99_proteome", count_fasta(snakemake.input.nr99)))

for track in ["SP_pos", "SP_neg"]:
    prefix = f"Secretome ({'+' if track == 'SP_pos' else '-'})"
    rows.append((f"{prefix}_after_SP_split",   count_fasta(f"{base}/02_sp_split/{sample}_{track}.faa")))
    rows.append((f"{prefix}_after_DeepLoc",    count_fasta(f"{base}/03_deeploc_filter/{sample}_{track}_extracellular.faa")))
    rows.append((f"{prefix}_after_TM_filter",  count_fasta(f"{base}/04_tm_filter/{sample}_{track}_noTM.faa")))

    if organism_type != "prokaryote":
        tp_path = f"{base}/05_targetp_filter/{sample}_{track}_noTP.faa"
        if os.path.exists(tp_path):
            rows.append((f"{prefix}_after_TP_filter", count_fasta(tp_path)))

        asafind_path = f"{base}/06_asafind_filter/{sample}_{track}_noPlastid.faa"
        if os.path.exists(asafind_path):
            rows.append((f"{prefix}_after_ASAFind_filter", count_fasta(asafind_path)))

    rows.append((f"{prefix}_final_secretome",  count_fasta(f"{base}/07_secretome/{sample}_{track}_secretome.faa")))

with open(snakemake.output.tsv, "w") as fh:
    fh.write("step\tcount\n")
    for step, count in rows:
        fh.write(f"{step}\t{count}\n")

print(f"Stats written to {snakemake.output.tsv}", file=sys.stderr)
