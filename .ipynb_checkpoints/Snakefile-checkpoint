# ==============================================================================
# Secretome Prediction Pipeline (Snakemake)
# ==============================================================================
# Linear filtering pipeline:
#   Step 1: CD-HIT         (redundancy reduction)       [secretome_pipeline env]
#   Step 2: SignalP 6.0    (split → SP+ / SP-)          [signalp6_env]
#   Step 3: DeepLoc 2.1    (retain extracellular)       [deeploc2 env]
#   Step 4: DeepTMHMM      (exclude transmembrane)      [apptainer / local / biolib]
#   Step 5: TargetP 2.0    (exclude organelle-targeted) [targetp2_env]       (eukaryotes only)
#   Step 6: ASAFind 2.0    (exclude complex-plastid)    [targetp2_env + venv] (plant_complex only)
#   Step 7: ER-retention   (exclude ER-retained)        [secretome_pipeline env] (eukaryotes only)
#
# Organism strategies:
#   prokaryote    — Steps 1–4 only (skip TargetP, ASAFind, ER)
#   animal        — skip Steps 5 cTP/luTP and Step 6
#   plant_simple  — TargetP cTP removal, skip Step 6
#   plant_complex — TargetP cTP/luTP + ASAFind
#
# Usage:
#   conda activate secretome_pipeline
#   snakemake --cores 8 --configfile config.yaml
#   snakemake --cores 8 --configfile config.yaml -np   # dry run
# ==============================================================================

import os

configfile: "config.yaml"

SAMPLES = list(config["samples"].keys())
TRACKS  = ["SP_pos", "SP_neg"]

# ---------------------------------------------------------------------------
# Conda activation helper
# ---------------------------------------------------------------------------
CONDA_BASE = config["tools"]["conda_base"]
CONDA_ACT  = f"source {CONDA_BASE}/etc/profile.d/conda.sh && conda activate"


def get_organism_type(sample):
    return config["samples"][sample]["organism_type"]


def is_prokaryote(sample):
    return get_organism_type(sample) == "prokaryote"


def get_targetp_classes(sample):
    otype = get_organism_type(sample)
    return config["targetp_remove"].get(otype, [])


def needs_asafind(sample):
    return get_organism_type(sample) == "plant_complex"


def get_targetp_org(sample):
    """TargetP organism flag: 'pl' for simple-plastid plants, 'non-pl' otherwise.
    Complex-plastid organisms use non-pl because their plastid targeting is
    detected by ASAFind, not by TargetP's plant model.
    Not called for prokaryotes (TargetP is skipped)."""
    otype = get_organism_type(sample)
    if otype == "plant_simple":
        return "pl"
    return "non-pl"


# ---------------------------------------------------------------------------
# Input function: resolve the file that feeds into the secretome output.
# Prokaryotes:    output of TM filtering (Step 4) — skip TargetP/ASAFind/ER
# plant_complex:  output of ASAFind filtering (Step 6) → ER filter
# Others:         output of TargetP filtering (Step 5) → ER filter
# ---------------------------------------------------------------------------
def er_input_faa(wildcards):
    if is_prokaryote(wildcards.sample):
        return f"results/{wildcards.sample}/04_tm_filter/{wildcards.sample}_{wildcards.track}_noTM.faa"
    elif needs_asafind(wildcards.sample):
        return f"results/{wildcards.sample}/06_asafind_filter/{wildcards.sample}_{wildcards.track}_noPlastid.faa"
    else:
        return f"results/{wildcards.sample}/05_targetp_filter/{wildcards.sample}_{wildcards.track}_noTP.faa"


def er_input_csv(wildcards):
    if is_prokaryote(wildcards.sample):
        return f"results/{wildcards.sample}/04_tm_filter/{wildcards.sample}_{wildcards.track}_noTM.csv"
    elif needs_asafind(wildcards.sample):
        return f"results/{wildcards.sample}/06_asafind_filter/{wildcards.sample}_{wildcards.track}_noPlastid.csv"
    else:
        return f"results/{wildcards.sample}/05_targetp_filter/{wildcards.sample}_{wildcards.track}_noTP.csv"


# ===========================================================================
# Target rule
# ===========================================================================
rule all:
    input:
        expand(
            "results/{sample}/07_secretome/{sample}_{track}_secretome.faa",
            sample=SAMPLES,
            track=TRACKS,
        ),
        expand(
            "results/{sample}/07_secretome/{sample}_{track}_secretome.csv",
            sample=SAMPLES,
            track=TRACKS,
        ),
        expand(
            "results/{sample}/summary/{sample}_pipeline_stats.tsv",
            sample=SAMPLES,
        ),


# ===========================================================================
# Step 1 — CD-HIT redundancy reduction  [secretome_pipeline env]
# ===========================================================================
rule cdhit:
    input:
        faa=lambda wc: config["samples"][wc.sample]["proteome"],
    output:
        faa="results/{sample}/01_cdhit/{sample}_nr99.faa",
        clstr="results/{sample}/01_cdhit/{sample}_nr99.faa.clstr",
    params:
        identity=config["cdhit"]["identity"],
        word_size=config["cdhit"]["word_size"],
    threads: config["cdhit"]["threads"]
    log:
        "logs/{sample}/01_cdhit.log",
    shell:
        """
        cd-hit -i {input.faa} -o {output.faa} \
            -c {params.identity} -n {params.word_size} \
            -T {threads} -M 0 \
            > {log} 2>&1
        """


# ===========================================================================
# Step 2a — SignalP 6.0  [signalp6_env]
# ===========================================================================
rule signalp:
    input:
        faa="results/{sample}/01_cdhit/{sample}_nr99.faa",
    output:
        txt="results/{sample}/02_signalp/{sample}_signalp_summary.signalp5",
        outdir=directory("results/{sample}/02_signalp/output"),
    params:
        organism=lambda wc: config["samples"][wc.sample]["signalp_organism"],
        conda_env=config["tools"]["signalp"]["conda_env"],
        cmd=config["tools"]["signalp"]["cmd"],
    threads: config["threads"]
    log:
        "logs/{sample}/02_signalp.log",
    shell:
        """
        {CONDA_ACT} {params.conda_env}
        {params.cmd} --fastafile {input.faa} \
            --organism {params.organism} \
            --output_dir {output.outdir} \
            --format txt \
            --mode fast \
            > {log} 2>&1
        cp {output.outdir}/prediction_results.txt {output.txt}
        """


# ===========================================================================
# Step 2b — Split into SP+ / SP-  [secretome_pipeline env]
# ===========================================================================
rule split_by_sp:
    input:
        faa="results/{sample}/01_cdhit/{sample}_nr99.faa",
        signalp="results/{sample}/02_signalp/{sample}_signalp_summary.signalp5",
    output:
        sp_pos_faa="results/{sample}/02_sp_split/{sample}_SP_pos.faa",
        sp_neg_faa="results/{sample}/02_sp_split/{sample}_SP_neg.faa",
        sp_pos_csv="results/{sample}/02_sp_split/{sample}_SP_pos.csv",
        sp_neg_csv="results/{sample}/02_sp_split/{sample}_SP_neg.csv",
    log:
        "logs/{sample}/02_split_by_sp.log",
    script:
        "workflow/scripts/split_by_sp.py"


# ===========================================================================
# Step 3a — DeepLoc 2.1  [deeploc2 env]
# ===========================================================================
rule deeploc:
    input:
        faa="results/{sample}/02_sp_split/{sample}_{track}.faa",
    output:
        csv="results/{sample}/03_deeploc/{sample}_{track}_deeploc.csv",
    params:
        conda_env=config["tools"]["deeploc"]["conda_env"],
        cmd=config["tools"]["deeploc"]["cmd"],
        model=config["tools"]["deeploc"]["model"],
        device=config["tools"]["deeploc"]["device"],
        outdir="results/{sample}/03_deeploc/{sample}_{track}_deeploc_out",
    log:
        "logs/{sample}/03_deeploc_{track}.log",
    resources:
        gpu=1,
    shell:
        """
        {CONDA_ACT} {params.conda_env}
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        {params.cmd} --fasta {input.faa} \
            --output {params.outdir} \
            --model {params.model} \
            --device {params.device} \
            > {log} 2>&1
        cp {params.outdir}/*.csv {output.csv}
        """


# ===========================================================================
# Step 3b — Retain extracellular proteins  [secretome_pipeline env]
# ===========================================================================
rule filter_deeploc:
    input:
        faa="results/{sample}/02_sp_split/{sample}_{track}.faa",
        csv="results/{sample}/02_sp_split/{sample}_{track}.csv",
        deeploc="results/{sample}/03_deeploc/{sample}_{track}_deeploc.csv",
    output:
        faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        csv="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.csv",
    log:
        "logs/{sample}/03_filter_deeploc_{track}.log",
    script:
        "workflow/scripts/filter_deeploc.py"


# ===========================================================================
# Step 4a — DeepTMHMM  [apptainer / local / biolib_local / biolib remote]
# ===========================================================================
DEEPTMHMM_METHOD = config["tools"]["deeptmhmm"]["method"]

if DEEPTMHMM_METHOD == "apptainer":
    rule deeptmhmm:
        input:
            faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        output:
            gff3="results/{sample}/04_deeptmhmm/{sample}_{track}_TMRs.gff3",
            outdir=directory("results/{sample}/04_deeptmhmm/{track}_output"),
        params:
            sif=config["tools"]["deeptmhmm"]["sif"],
        log:
            "logs/{sample}/04_deeptmhmm_{track}.log",
        shell:
            """
            mkdir -p {output.outdir}
            apptainer exec --nv --bind "$PWD:/work" {params.sif} \
                bash -lc "cd /work/{output.outdir} && \
                python3 /openprotein/predict.py --fasta /work/{input.faa}" \
                > {log} 2>&1
            cp {output.outdir}/TMRs.gff3 {output.gff3}
            """
elif DEEPTMHMM_METHOD == "local":
    rule deeptmhmm:
        input:
            faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        output:
            gff3="results/{sample}/04_deeptmhmm/{sample}_{track}_TMRs.gff3",
            outdir=directory("results/{sample}/04_deeptmhmm/{track}_output"),
        params:
            conda_env=config["tools"]["deeptmhmm"]["conda_env"],
            predict_py=config["tools"]["deeptmhmm"]["predict_py"],
        log:
            "logs/{sample}/04_deeptmhmm_{track}.log",
        shell:
            """
            WORKDIR=$(pwd)
            {CONDA_ACT} {params.conda_env}
            cd $(dirname {params.predict_py})
            python predict.py \
                --fasta $WORKDIR/{input.faa} \
                --output-dir $WORKDIR/{output.outdir} \
                > $WORKDIR/{log} 2>&1
            cp $WORKDIR/{output.outdir}/TMRs.gff3 $WORKDIR/{output.gff3}
            """
elif DEEPTMHMM_METHOD == "biolib_local":
    rule deeptmhmm:
        input:
            faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        output:
            gff3="results/{sample}/04_deeptmhmm/{sample}_{track}_TMRs.gff3",
            outdir=directory("results/{sample}/04_deeptmhmm/{track}_output"),
        log:
            "logs/{sample}/04_deeptmhmm_{track}.log",
        shell:
            """
            mkdir -p {output.outdir}
            FASTA_ABS=$(cd $(dirname {input.faa}) && pwd)/$(basename {input.faa})
            cd {output.outdir}
            biolib run --local 'DTU/DeepTMHMM:1.0.24' \
                --fasta "$FASTA_ABS" \
                > ../../$(dirname {log})/$(basename {log}) 2>&1
            cp TMRs.gff3 ../$(basename {output.gff3}) || \
                cp biolib_results/TMRs.gff3 ../$(basename {output.gff3})
            """
else:
    rule deeptmhmm:
        input:
            faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        output:
            gff3="results/{sample}/04_deeptmhmm/{sample}_{track}_TMRs.gff3",
            outdir=directory("results/{sample}/04_deeptmhmm/{track}_output"),
        params:
            conda_env=config["tools"]["deeptmhmm"].get("conda_env", "deeptmhmm_env"),
        log:
            "logs/{sample}/04_deeptmhmm_{track}.log",
        shell:
            """
            {CONDA_ACT} {params.conda_env}
            biolib run DTU/DeepTMHMM --fasta {input.faa} \
                --output {output.outdir} \
                > {log} 2>&1
            cp {output.outdir}/TMRs.gff3 {output.gff3}
            """


# ===========================================================================
# Step 4b — Filter TM proteins  [secretome_pipeline env]
# ===========================================================================
rule filter_tm:
    input:
        faa="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.faa",
        csv="results/{sample}/03_deeploc_filter/{sample}_{track}_extracellular.csv",
        gff3="results/{sample}/04_deeptmhmm/{sample}_{track}_TMRs.gff3",
    output:
        faa="results/{sample}/04_tm_filter/{sample}_{track}_noTM.faa",
        csv="results/{sample}/04_tm_filter/{sample}_{track}_noTM.csv",
        excluded="results/{sample}/04_tm_filter/{sample}_{track}_TM_excluded.txt",
    log:
        "logs/{sample}/04_filter_tm_{track}.log",
    script:
        "workflow/scripts/filter_tm.py"


# ===========================================================================
# Step 5a — TargetP 2.0  [targetp2_env]
# ===========================================================================
rule targetp:
    input:
        faa="results/{sample}/04_tm_filter/{sample}_{track}_noTM.faa",
    output:
        txt="results/{sample}/05_targetp/{sample}_{track}_targetp.txt",
        outdir=directory("results/{sample}/05_targetp/{track}_output"),
    params:
        conda_env=config["tools"]["targetp"]["conda_env"],
        cmd=config["tools"]["targetp"]["cmd"],
        org=lambda wc: get_targetp_org(wc.sample),
    log:
        "logs/{sample}/05_targetp_{track}.log",
    shell:
        """
        {CONDA_ACT} {params.conda_env}
        mkdir -p {output.outdir}
        {params.cmd} -fasta {input.faa} -org {params.org} \
            -prefix {output.outdir}/{wildcards.sample}_{wildcards.track} \
            > {log} 2>&1
        cp {output.outdir}/{wildcards.sample}_{wildcards.track}_summary.targetp2 {output.txt}
        """


# ===========================================================================
# Step 5b — Filter TargetP predictions  [secretome_pipeline env]
# ===========================================================================
rule filter_targetp:
    input:
        faa="results/{sample}/04_tm_filter/{sample}_{track}_noTM.faa",
        csv="results/{sample}/04_tm_filter/{sample}_{track}_noTM.csv",
        targetp="results/{sample}/05_targetp/{sample}_{track}_targetp.txt",
    output:
        faa="results/{sample}/05_targetp_filter/{sample}_{track}_noTP.faa",
        csv="results/{sample}/05_targetp_filter/{sample}_{track}_noTP.csv",
        excluded="results/{sample}/05_targetp_filter/{sample}_{track}_TP_excluded.txt",
    params:
        remove_classes=lambda wc: get_targetp_classes(wc.sample),
    log:
        "logs/{sample}/05_filter_targetp_{track}.log",
    script:
        "workflow/scripts/filter_targetp.py"


# ===========================================================================
# Step 6a — ASAFind 2.0 (complex plastid only)  [targetp2_env + venv]
#           Requires TargetP output (non-pl model) to detect bipartite
#           plastid-targeting signals in complex-plastid organisms.
# ===========================================================================
rule asafind:
    input:
        faa="results/{sample}/05_targetp_filter/{sample}_{track}_noTP.faa",
        signalp="results/{sample}/02_signalp/{sample}_signalp_summary.signalp5",
        targetp="results/{sample}/05_targetp/{sample}_{track}_targetp.txt",
    output:
        txt="results/{sample}/06_asafind/{sample}_{track}_asafind.txt",
    params:
        conda_env=config["tools"]["asafind"]["conda_env"],
        venv=config["tools"]["asafind"]["venv"],
        script=config["tools"]["asafind"]["script"],
        asafind_dir=os.path.dirname(config["tools"]["asafind"]["script"]),
    log:
        "logs/{sample}/06_asafind_{track}.log",
    resources:
        asafind_lock=1,
    shell:
        """
        WORKDIR=$(pwd)
        {CONDA_ACT} {params.conda_env}
        source {params.venv}/bin/activate
        cd {params.asafind_dir}
        rm -rf temp/ output/
        python S0_ASAFind.py \
            -f $WORKDIR/{input.faa} \
            -p $WORKDIR/{input.targetp} \
            > $WORKDIR/{log} 2>&1
        mkdir -p $WORKDIR/$(dirname {output.txt})
        unzip -o output/output_*.zip "*.tab" -d output/ >> $WORKDIR/{log} 2>&1
        cp output/*.tab $WORKDIR/{output.txt}
        """


# ===========================================================================
# Step 6b — Filter ASAFind plastid-targeted  [secretome_pipeline env]
# ===========================================================================
rule filter_asafind:
    input:
        faa="results/{sample}/05_targetp_filter/{sample}_{track}_noTP.faa",
        csv="results/{sample}/05_targetp_filter/{sample}_{track}_noTP.csv",
        asafind="results/{sample}/06_asafind/{sample}_{track}_asafind.txt",
    output:
        faa="results/{sample}/06_asafind_filter/{sample}_{track}_noPlastid.faa",
        csv="results/{sample}/06_asafind_filter/{sample}_{track}_noPlastid.csv",
        excluded="results/{sample}/06_asafind_filter/{sample}_{track}_plastid_excluded.txt",
    log:
        "logs/{sample}/06_filter_asafind_{track}.log",
    script:
        "workflow/scripts/filter_asafind.py"


# ===========================================================================
# Step 7 — ER-retention motif exclusion → FINAL SECRETOME
#          [secretome_pipeline env]
# ===========================================================================
rule filter_er:
    input:
        faa=er_input_faa,
        csv=er_input_csv,
    output:
        faa="results/{sample}/07_secretome/{sample}_{track}_secretome.faa",
        csv="results/{sample}/07_secretome/{sample}_{track}_secretome.csv",
        excluded="results/{sample}/07_secretome/{sample}_{track}_ER_excluded.txt",
    params:
        motifs=config["er_motifs"],
    log:
        "logs/{sample}/07_filter_er_{track}.log",
    script:
        "workflow/scripts/filter_er.py"


# ===========================================================================
# Summary statistics  [secretome_pipeline env]
# ===========================================================================
rule pipeline_stats:
    input:
        nr99="results/{sample}/01_cdhit/{sample}_nr99.faa",
        sp_pos="results/{sample}/07_secretome/{sample}_SP_pos_secretome.faa",
        sp_neg="results/{sample}/07_secretome/{sample}_SP_neg_secretome.faa",
    output:
        tsv="results/{sample}/summary/{sample}_pipeline_stats.tsv",
    log:
        "logs/{sample}/pipeline_stats.log",
    script:
        "workflow/scripts/pipeline_stats.py"
