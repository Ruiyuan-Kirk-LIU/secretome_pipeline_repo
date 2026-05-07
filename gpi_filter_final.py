#!/usr/bin/env python3
"""
Optional post-pipeline step: GPI-anchor filtering for Secretome (SP+).

This script filters GPI-anchored proteins from the SP+ secretome using
NetGPI web server output, and organises final secretomes into 09_final/.

Usage:
    # Process a single sample:
    python gpi_filter_final.py --sample Smic

    # Process multiple samples:
    python gpi_filter_final.py --sample Smic Cgor Dtre

    # Process all samples found in results/:
    python gpi_filter_final.py --all

Prerequisites:
    1. Run the main Snakemake pipeline to completion (07_secretome/ must exist).
    2. Submit each sample's SP+ FASTA from 02_sp_split/ to NetGPI 1.1:
           https://services.healthtech.dtu.dk/services/NetGPI-1.1/
    3. Place each output file in:
           results/{sample}/08_gpi/{sample}_gpi_output.txt
       The expected filename pattern is {sample}_gpi_output*.txt

Directory structure after running this script:
    results/{sample}/
    ├── 07_secretome/
    │   ├── {sample}_SP_pos_secretome.faa     (from pipeline)
    │   └── {sample}_SP_neg_secretome.faa     (from pipeline)
    ├── 08_gpi_filter/
    │   ├── {sample}_SP_pos_noGPI.faa         (SP+ minus GPI-anchored)
    │   ├── {sample}_SP_pos_noGPI.csv
    │   ├── {sample}_SP_pos_GPI_excluded.faa  (GPI-anchored proteins)
    │   └── {sample}_SP_pos_GPI_excluded.txt  (excluded IDs)
    └── 09_final/
        ├── {sample}_SP_pos_secretome_final.faa
        ├── {sample}_SP_pos_secretome_final.csv
        ├── {sample}_SP_neg_secretome_final.faa
        ├── {sample}_SP_neg_secretome_final.csv
        └── {sample}_pipeline_stats_final.tsv
"""

import argparse
import glob
import os
import shutil
import sys

import pandas as pd
from Bio import SeqIO


def find_gpi_file(sample, results_dir="results"):
    """Find the NetGPI output file for a sample."""
    gpi_dir = os.path.join(results_dir, sample, "08_gpi")
    if not os.path.isdir(gpi_dir):
        return None
    candidates = glob.glob(os.path.join(gpi_dir, f"{sample}_gpi_output*.txt"))
    if not candidates:
        # Try any .txt file in the directory
        candidates = glob.glob(os.path.join(gpi_dir, "*.txt"))
    if len(candidates) == 1:
        return candidates[0]
    elif len(candidates) > 1:
        print(f"  [WARNING] Multiple GPI files found for {sample}, using: {candidates[0]}")
        return candidates[0]
    return None


def parse_netgpi(gpi_file):
    """Parse NetGPI output and return sets of GPI-anchored and non-anchored IDs."""
    gpi_ids = set()
    non_gpi_ids = set()
    with open(gpi_file) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            pid = parts[0].split()[0]
            pred = parts[2].strip()
            if pred == "GPI-Anchored":
                gpi_ids.add(pid)
            else:
                non_gpi_ids.add(pid)
    return gpi_ids, non_gpi_ids


def process_sample(sample, results_dir="results"):
    """Run GPI filtering and produce final secretomes for one sample."""
    base = os.path.join(results_dir, sample)

    # --- Check prerequisites ---
    sp_pos_faa = os.path.join(base, "07_secretome", f"{sample}_SP_pos_secretome.faa")
    sp_pos_csv = os.path.join(base, "07_secretome", f"{sample}_SP_pos_secretome.csv")
    sp_neg_faa = os.path.join(base, "07_secretome", f"{sample}_SP_neg_secretome.faa")
    sp_neg_csv = os.path.join(base, "07_secretome", f"{sample}_SP_neg_secretome.csv")

    if not os.path.exists(sp_pos_faa):
        print(f"  [SKIP] {sample}: SP+ secretome not found at {sp_pos_faa}")
        return False

    gpi_file = find_gpi_file(sample, results_dir)
    if gpi_file is None:
        print(f"  [SKIP] {sample}: No NetGPI output found in {base}/08_gpi/")
        print(f"         Place the file as: {base}/08_gpi/{sample}_gpi_output.txt")
        return False

    print(f"  NetGPI file: {gpi_file}")

    # --- Parse NetGPI ---
    gpi_ids, non_gpi_ids = parse_netgpi(gpi_file)
    print(f"  NetGPI: {len(gpi_ids)} GPI-anchored, {len(non_gpi_ids)} not anchored")

    # --- Filter SP+ secretome ---
    gpi_filter_dir = os.path.join(base, "08_gpi_filter")
    os.makedirs(gpi_filter_dir, exist_ok=True)

    records = list(SeqIO.parse(sp_pos_faa, "fasta"))
    kept = []
    excluded = []
    for rec in records:
        rid = rec.id.split()[0]
        if rid in gpi_ids:
            excluded.append(rec)
        else:
            kept.append(rec)

    # Write filtered FASTA
    noGPI_faa = os.path.join(gpi_filter_dir, f"{sample}_SP_pos_noGPI.faa")
    SeqIO.write(kept, noGPI_faa, "fasta")

    # Write excluded FASTA
    excluded_faa = os.path.join(gpi_filter_dir, f"{sample}_SP_pos_GPI_excluded.faa")
    SeqIO.write(excluded, excluded_faa, "fasta")

    # Write excluded IDs
    excluded_txt = os.path.join(gpi_filter_dir, f"{sample}_SP_pos_GPI_excluded.txt")
    with open(excluded_txt, "w") as fh:
        for rec in sorted(excluded, key=lambda r: r.id):
            fh.write(rec.id.split()[0] + "\n")

    # Write filtered CSV
    noGPI_csv = os.path.join(gpi_filter_dir, f"{sample}_SP_pos_noGPI.csv")
    if os.path.exists(sp_pos_csv):
        df = pd.read_csv(sp_pos_csv)
        df["Protein_ID"] = df["Protein_ID"].astype(str)
        df_kept = df[~df["Protein_ID"].isin(gpi_ids)]
        df_kept.to_csv(noGPI_csv, index=False)
    else:
        pd.DataFrame({"Protein_ID": [r.id.split()[0] for r in kept]}).to_csv(noGPI_csv, index=False)

    print(f"  SP+ secretome: {len(records)} → {len(kept)} kept, {len(excluded)} GPI removed")
    print(f"  → {noGPI_faa}")
    print(f"  → {excluded_faa}")

    # --- Produce 09_final ---
    final_dir = os.path.join(base, "09_final")
    os.makedirs(final_dir, exist_ok=True)

    # SP+ final = GPI-filtered
    sp_pos_final_faa = os.path.join(final_dir, f"{sample}_SP_pos_secretome_final.faa")
    sp_pos_final_csv = os.path.join(final_dir, f"{sample}_SP_pos_secretome_final.csv")
    shutil.copy2(noGPI_faa, sp_pos_final_faa)
    shutil.copy2(noGPI_csv, sp_pos_final_csv)

    # SP- final = copy from 07_secretome
    sp_neg_final_faa = os.path.join(final_dir, f"{sample}_SP_neg_secretome_final.faa")
    sp_neg_final_csv = os.path.join(final_dir, f"{sample}_SP_neg_secretome_final.csv")
    if os.path.exists(sp_neg_faa):
        shutil.copy2(sp_neg_faa, sp_neg_final_faa)
    else:
        open(sp_neg_final_faa, "w").close()
    if os.path.exists(sp_neg_csv):
        shutil.copy2(sp_neg_csv, sp_neg_final_csv)
    else:
        pd.DataFrame(columns=["Protein_ID"]).to_csv(sp_neg_final_csv, index=False)

    print(f"  → {sp_pos_final_faa}")
    print(f"  → {sp_neg_final_faa}")

    # --- Update summary stats ---
    stats_file = os.path.join(base, "summary", f"{sample}_pipeline_stats.tsv")
    final_stats = os.path.join(final_dir, f"{sample}_pipeline_stats_final.tsv")

    if os.path.exists(stats_file):
        stats = pd.read_csv(stats_file, sep="\t")
    else:
        stats = pd.DataFrame(columns=["step", "count"])

    # Count final sequences
    sp_pos_final_count = sum(1 for _ in SeqIO.parse(sp_pos_final_faa, "fasta")) if os.path.getsize(sp_pos_final_faa) > 0 else 0
    sp_neg_final_count = sum(1 for _ in SeqIO.parse(sp_neg_final_faa, "fasta")) if os.path.getsize(sp_neg_final_faa) > 0 else 0

    new_rows = pd.DataFrame([
        {"step": "Secretome (+)_after_GPI_filter", "count": len(kept)},
        {"step": "Secretome (+)_final_with_GPI", "count": sp_pos_final_count},
        {"step": "Secretome (-)_final_with_GPI", "count": sp_neg_final_count},
    ])
    stats_final = pd.concat([stats, new_rows], ignore_index=True)
    stats_final.to_csv(final_stats, sep="\t", index=False)

    print(f"  → {final_stats}")
    print(f"  SP+ final: {sp_pos_final_count} | SP- final: {sp_neg_final_count}")
    return True


def find_all_samples(results_dir="results"):
    """Find all samples that have completed the pipeline."""
    samples = []
    if os.path.isdir(results_dir):
        for d in sorted(os.listdir(results_dir)):
            secretome_dir = os.path.join(results_dir, d, "07_secretome")
            if os.path.isdir(secretome_dir):
                samples.append(d)
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Post-pipeline GPI filtering and final secretome generation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gpi_filter_final.py --sample Smic
  python gpi_filter_final.py --sample Smic Cgor Dtre
  python gpi_filter_final.py --all
  python gpi_filter_final.py --all --results-dir /path/to/results
        """,
    )
    parser.add_argument("--sample", nargs="+", help="Sample name(s) to process")
    parser.add_argument("--all", action="store_true", help="Process all samples in results/")
    parser.add_argument("--results-dir", default="results", help="Path to results directory (default: results)")
    args = parser.parse_args()

    if not args.sample and not args.all:
        parser.print_help()
        sys.exit(1)

    if args.all:
        samples = find_all_samples(args.results_dir)
        if not samples:
            print(f"No completed samples found in {args.results_dir}/")
            sys.exit(1)
        print(f"Found {len(samples)} samples: {', '.join(samples)}")
    else:
        samples = args.sample

    success = 0
    skipped = 0
    for sample in samples:
        print(f"\n===== {sample} =====")
        if process_sample(sample, args.results_dir):
            success += 1
        else:
            skipped += 1

    print(f"\n{'='*40}")
    print(f"Processed: {success} | Skipped: {skipped}")
    if skipped > 0:
        print("Skipped samples need NetGPI output in results/{sample}/08_gpi/")


if __name__ == "__main__":
    main()
