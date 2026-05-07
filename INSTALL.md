# Installation Guide

This pipeline requires five conda environments plus an ASAFind Python venv. Most tools are licensed from DTU and must be downloaded separately after accepting their academic license.

## Prerequisites

- Conda (Miniforge/Miniconda/Anaconda)
- Git
- For HPC: PBS scheduler, Apptainer (optional for DeepTMHMM)
- For Mac: Docker Desktop (optional for DeepTMHMM via biolib)

## Environment Setup

### Environment 0: Pipeline Runner

Contains Snakemake, CD-HIT, and Python dependencies for the filtering scripts.

```bash
conda create -n secretome_pipeline -c conda-forge -c bioconda \
    snakemake cd-hit --strict-channel-priority
conda activate secretome_pipeline
pip install biopython pandas
```

### Environment 1: SignalP 6.0

Requires a license from DTU. Download from: https://services.healthtech.dtu.dk/services/SignalP-6.0/

```bash
conda create -n signalp6_env python=3.9 "numpy<2"
conda activate signalp6_env

# Install from the downloaded package
cd signalp6_fast/signalp-6-package
pip install .

# Copy model weights to the installed location
SIGNALP_DIR=$(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))")
cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/

# Test
signalp6 --fastafile test.fasta --organism other \
    --output_dir test_out --format txt --mode fast
```

**Known issue**: `numpy>=2` causes `RuntimeError: Numpy is not available`. Always pin `numpy<2`.

### Environment 2: DeepLoc 2.1

Requires a license from DTU. Download from: https://services.healthtech.dtu.dk/services/DeepLoc-2.1/

```bash
conda create -y -n deeploc2 -c conda-forge \
    python=3.10 "numpy<2" onnxruntime=1.17 pip
conda activate deeploc2
pip install deeploc-2.1.All.tar.gz

# Test
deeploc2 --fasta test.fasta --output test_out --model Fast --device cpu
```

**Notes**:
- `--output` takes a **directory**, not a file path
- Use `--device cuda` on GPU nodes, `--device mps` on Apple Silicon, `--device cpu` elsewhere
- GPU jobs should be limited to one at a time (`--resources gpu=1`) to prevent OOM

### Environment 3: TargetP 2.0 & ASAFind 2.0

#### TargetP 2.0

Download the platform-specific package from DTU: https://services.healthtech.dtu.dk/services/TargetP-2.0/

- Linux: `targetp-2.0.Linux.tar.gz`
- macOS: `targetp-2.0.Darwin.tar.gz`

```bash
conda create -n targetp2_env python=3.9 "numpy<2"
conda activate targetp2_env

# Extract
tar -xvzf targetp-2.0.tar.gz

# Test
targetp-2.0/bin/targetp -fasta test.fasta -org non-pl -format short
```

**macOS Gatekeeper**: If blocked, run:
```bash
sudo xattr -cr targetp-2.0/
```

#### ASAFind 2.0

```bash
# Clone the repository
git clone https://github.com/ASAFind/ASAFind-2.git

# Create a Python venv (inside the ASAFind directory or alongside it)
cd ASAFind-2   # or your preferred location
python3 -m venv asafind_command_line
source asafind_command_line/bin/activate
pip install pandas

# Test (must run from the ASAFind-2 directory)
python S0_ASAFind.py -f test.fasta -p test_summary.targetp2
```

**Notes**:
- ASAFind must be run from its own directory (uses relative imports for `S1_ASAFind.py`, `fill_constants.py`)
- Output is a ZIP file in `output/` containing a `.tab` file
- The pipeline handles this via `cd` and `unzip` in the shell command
- Concurrent runs will clobber `temp/` and `output/` — the pipeline uses `asafind_lock=1` to serialize

### Environment 4: DeepTMHMM

Multiple installation options:

#### Option A: DTU Academic Package (recommended for local)

Download from: https://services.healthtech.dtu.dk/services/DeepTMHMM-1.0/

```bash
conda create -n deeptmhmm_env python=3.8 -y
conda activate deeptmhmm_env
pip install wheel Cython==0.29.37 pkgconfig==1.5.5

# Install PyTorch (CPU version; use CUDA wheel for GPU)
pip install torch

# Install dependencies
cd DeepTMHMM-Academic-License-v1.0
pip install -r requirements.txt

# Test
python predict.py --fasta sample.fasta --output-dir result1
```

Set `method: "local"` in config.yaml.

#### Option B: Apptainer/Singularity (recommended for HPC)

Build or obtain the DeepTMHMM container image.

```bash
# Test
apptainer exec --nv --bind "$PWD:/work" deeptmhmm \
    bash -lc "cd /work && python3 /openprotein/predict.py --fasta test.fasta"
```

Set `method: "apptainer"` in config.yaml.

#### Option C: biolib (requires Docker or internet)

```bash
conda create -n deeptmhmm_env python=3.12 -y
conda activate deeptmhmm_env
pip install pybiolib

# Local (requires Docker Desktop)
biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta test.fasta

# Remote (requires internet)
biolib run DTU/DeepTMHMM --fasta test.fasta
```

Set `method: "biolib_local"` or `method: "biolib"` in config.yaml.

## Verifying Installation

After installing all environments, verify Python dependencies for the filtering scripts:

```bash
conda activate secretome_pipeline
python -c "from Bio import SeqIO; import pandas; print('OK')"
```

Then run a dry test:

```bash
snakemake --cores 1 -np --configfile config.yaml --resources gpu=1 asafind_lock=1
```

## Troubleshooting

| Issue | Solution |
|---|---|
| `signalp6: command not found` | Conda activation failed inside shell rule — check `conda_base` path in config.yaml |
| `NumPy 2.0` errors with SignalP | Pin `numpy<2` in the signalp6_env |
| DeepLoc creates a directory instead of file | Expected — `--output` takes a directory, pipeline handles the `cp` |
| DeepLoc OOM / Killed | Run with `--resources gpu=1` to serialize GPU jobs; increase PBS memory |
| ASAFind `ModuleNotFoundError: fill_constants` | Must run from ASAFind's own directory — pipeline does this via `cd` |
| ASAFind `Directory not empty` race condition | Clean `temp/` and `output/` in ASAFind dir; ensure `asafind_lock=1` is passed |
| TargetP blocked on macOS | Run `sudo xattr -cr targetp-2.0/` |
| `nano` segfault on HPC | Use `vi` or `sed` instead |
| Conda verification errors / ClobberErrors | Run `conda clean --packages --tarballs` and use `--strict-channel-priority` |
