# Secretome Prediction Pipeline

A Snakemake pipeline for predicting high-confidence secreted proteins from eukaryotic proteomes. Supports three organism strategies and produces two secretome classes per sample: **Secretome (SP+)** and **Secretome (SP-)**.

Tested on HPC (PBS scheduler, NVIDIA A100) and Apple M2 Mac Studio.

## Pipeline Overview

```
Raw proteome
  │
  ├─ Step 1: CD-HIT         → redundancy reduction (99% identity)
  ├─ Step 2: SignalP 6.0    → split into SP+ (signal peptide) / SP- (no SP)
  ├─ Step 3: DeepLoc 2.1    → retain extracellular proteins
  ├─ Step 4: DeepTMHMM      → exclude transmembrane proteins
  ├─ Step 5: TargetP 2.0    → exclude organelle-targeted proteins
  ├─ Step 6: ASAFind 2.0    → exclude complex-plastid targeted (plant_complex only)
  ├─ Step 7: ER retention   → exclude ER-retained proteins (KDEL/HDEL motifs)
  │
  ├─ Secretome (SP+)  ← results/{sample}/07_secretome/{sample}_SP_pos_secretome.faa
  └─ Secretome (SP-)  ← results/{sample}/07_secretome/{sample}_SP_neg_secretome.faa

Optional post-pipeline:
  ├─ Step 8: NetGPI 1.1     → exclude GPI-anchored (SP+ only, via web server)
  ├─ 09_final/              → final secretomes with GPI filtering
```

## Organism Strategies

| Strategy | TargetP removes | ASAFind | TargetP `-org` | Example organisms |
|---|---|---|---|---|
| `animal` | mTP | skipped | `non-pl` | Metazoans |
| `plant_simple` | mTP, cTP | skipped | `pl` | Land plants, green algae |
| `plant_complex` | mTP, cTP, luTP | yes | `non-pl` | Dinoflagellates, diatoms, cryptophytes |

## Quick Start

```bash
# 1. Install environments (see INSTALL.md for details)
# 2. Place proteomes in data/raw/
# 3. Edit config/config_hpc.yaml or config/config_mac.yaml
# 4. Run
conda activate secretome_pipeline
snakemake --cores 8 --configfile config/config_hpc.yaml \
    --resources gpu=1 asafind_lock=1
```

## Configuration

Copy and edit the appropriate template:

```bash
# HPC with PBS
cp config/config_hpc.yaml config.yaml

# Local Mac
cp config/config_mac.yaml config.yaml
```

Add samples to the `samples:` section:

```yaml
samples:
  Cgor:
    proteome: "data/raw/Cgor.faa"
    organism_type: "plant_complex"
    signalp_organism: "other"
```

### DeepTMHMM Methods

The pipeline supports four methods for running DeepTMHMM, configured via `tools.deeptmhmm.method`:

| Method | Config key | Use case |
|---|---|---|
| `apptainer` | `sif:` path to image | HPC without internet |
| `local` | `conda_env:`, `predict_py:` | DTU academic package (any system) |
| `biolib_local` | (none) | Docker via biolib (Mac/Linux) |
| `biolib` | `conda_env:` | Remote biolib (needs internet) |

## Running on HPC (PBS)

```bash
# PBS job script (run_pipeline.pbs)
#!/bin/bash
#PBS -N secretome
#PBS -j oe
#PBS -l select=1:ncpus=16:ngpus=1:mem=64gb
#PBS -l walltime=100:00:00
set -eo pipefail
cd "${PBS_O_WORKDIR:?}"
module load conda/20260225
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate secretome_pipeline
snakemake --cores 8 --configfile config.yaml \
    --resources gpu=1 asafind_lock=1

# Submit
qsub run_pipeline.pbs
```

### Resource notes

- `gpu=1` limits DeepLoc to one job at a time (prevents GPU OOM)
- `asafind_lock=1` serialises ASAFind runs (shares temp/output directories)
- Request `mem=64gb` or more for full proteomes; `mem=32gb` may OOM with multiple DeepLoc jobs

## Running Locally (Mac/Linux)

```bash
conda activate secretome_pipeline
snakemake --cores 8 --configfile config/config_mac.yaml \
    --resources gpu=1 asafind_lock=1
```

## Optional: GPI-Anchor Filtering (Post-Pipeline)

GPI filtering for SP+ proteins uses the NetGPI 1.1 web server since no standalone package is available.

```bash
# 1. Submit SP+ FASTAs to https://services.healthtech.dtu.dk/services/NetGPI-1.1/
# 2. Place results in results/{sample}/08_gpi/
mkdir -p results/Smic/08_gpi
cp downloaded_result.txt results/Smic/08_gpi/Smic_gpi_output.txt

# 3. Run the post-pipeline script
python gpi_filter_final.py --sample Smic
# Or process all samples:
python gpi_filter_final.py --all
```

This produces:

```
results/{sample}/
├── 08_gpi_filter/
│   ├── {sample}_SP_pos_noGPI.faa
│   ├── {sample}_SP_pos_GPI_excluded.faa
│   └── {sample}_SP_pos_GPI_excluded.txt
└── 09_final/
    ├── {sample}_SP_pos_secretome_final.faa
    ├── {sample}_SP_neg_secretome_final.faa
    └── {sample}_pipeline_stats_final.tsv
```

## Output Structure

```
results/{sample}/
├── 01_cdhit/               # Non-redundant proteome
├── 02_signalp/             # SignalP 6.0 results
├── 02_sp_split/            # SP+ and SP- FASTAs
├── 03_deeploc/             # DeepLoc raw output
├── 03_deeploc_filter/      # Extracellular proteins
├── 04_deeptmhmm/           # DeepTMHMM raw output
├── 04_tm_filter/           # After TM removal
├── 05_targetp/             # TargetP raw output
├── 05_targetp_filter/      # After TargetP removal
├── 06_asafind/             # ASAFind raw output (complex plastid only)
├── 06_asafind_filter/      # After plastid removal (complex plastid only)
├── 07_secretome/           # ★ Pipeline secretomes (SP+ and SP-)
├── 08_gpi/                 # NetGPI results (manual, optional)
├── 08_gpi_filter/          # After GPI removal (optional)
├── 09_final/               # Final secretomes with GPI (optional)
└── summary/                # Pipeline statistics
```

## Resuming After Failures

Snakemake automatically resumes from the last successful step:

```bash
# Just resubmit — completed steps are skipped
qsub run_pipeline.pbs
```

To rerun a specific step:

```bash
snakemake --cores 8 --configfile config.yaml \
    --resources gpu=1 asafind_lock=1 \
    --forcerun filter_er
```

## Adding New Samples

Add entries to `config.yaml` and rerun. Snakemake only processes new/changed samples:

```yaml
samples:
  NewSpecies:
    proteome: "data/raw/NewSpecies.faa"
    organism_type: "plant_complex"
    signalp_organism: "other"
```

## Citation

If you use this pipeline, please cite the tools used:

- **CD-HIT**: Fu et al. (2012) Bioinformatics 28:3150-3152
- **SignalP 6.0**: Teufel et al. (2022) Nature Biotechnology 40:1023-1025
- **DeepLoc 2.1**: Thumuluri et al. (2022) Nucleic Acids Research 50:W228-W234
- **DeepTMHMM**: Hallgren et al. (2022) bioRxiv doi:10.1101/2022.04.08.487609
- **TargetP 2.0**: Armenteros et al. (2019) Nature Plants 5:401-408
- **ASAFind 2.0**: Gruber et al. (2015) Plant Cell 27:2688-2698
- **NetGPI 1.1**: Gíslason et al. (2021) Bioinformatics 37:2836-2842

## License

MIT
