# Snakefile

configfile: "config.yaml"

import re

# Executables / environments (set in config.yaml; defaults assume they are on PATH)
RSCRIPT          = config.get("rscript", "Rscript")
ENV_SPAGCN       = config.get("env_spagcn", "env/SpaGCN")      # micromamba prefix for SpaGCN
ENV_STLEARN      = config.get("env_stlearn", "env/stlearn")    # micromamba prefix for stLearn
PY_SPATIALLEIDEN = config.get("python_spatialleiden", "python")
PYTHON           = config.get("python", "python")              # python with scanpy (prepare_anndata)

# ---------------------------------------------------------------------------
# Number of principal components (config.yaml: PC)
#   PC: default -> every method uses its own package default (PC_DEFAULTS)
#   PC: <int>   -> every method uses that number
# A method block may still set n_pcs explicitly (optimize_parameter.py does this
# when it tunes stLearn's n_pcs); that value wins for that method only.
# ---------------------------------------------------------------------------
PC_DEFAULTS = {
    "clustering":    10,  # louvain / louvain_refine / SLM / leiden: Seurat FindNeighbors(dims = 1:10)
    "BANKSY":        20,  # Banksy::runBanksyPCA(npcs = 20)
    "BayesSpace":    15,  # BayesSpace::spatialCluster(d = 15)
    "SpatialPCA":    20,  # SpatialPCA::SpatialPCA_EstimateLoading(SpatialPCnum = 20)
    "SpaGCN":        50,  # SpaGCN.SpaGCN.train(num_pcs = 50)
    "stlearn":       50,  # stlearn.em.run_pca(n_comps = 50)
    "SpatialLeiden": 50,  # scanpy.pp.pca default (50 components)
}
STEP3_DEFAULT_NPCS = 50   # Seurat RunPCA(npcs = 50): PCs computed in Step 3 (04_reduced.rds)

PC_SETTING = str(config.get("PC", "default")).strip()
if PC_SETTING.lower() != "default" and not (PC_SETTING.isdigit() and int(PC_SETTING) > 0):
    raise ValueError(f"config PC must be 'default' or a positive integer, got {PC_SETTING!r}")

def n_pcs(method):
    """Number of PCs used by `method` (a key of PC_DEFAULTS)."""
    block = config.get(method) or {}
    if block.get("n_pcs") is not None:
        return int(block["n_pcs"])
    if PC_SETTING.lower() == "default":
        return PC_DEFAULTS[method]
    return int(PC_SETTING)

def step3_npcs():
    """PCs computed in Step 3: enough for every method that reads them (Seurat clustering, BayesSpace)."""
    base = STEP3_DEFAULT_NPCS if PC_SETTING.lower() == "default" else int(PC_SETTING)
    return max(base, n_pcs("clustering"), n_pcs("BayesSpace"))

# Step 0 QC options passed to prepare_seurat.R and prepare_anndata.py (config.yaml: qc)
_QC = config.get("qc") or {}
QC_ARGS = (f"--remove_mt {str(bool(_QC.get('remove_mt', True))).upper()} "
           f"--mt_pattern '{_QC.get('mt_pattern', '^mt-')}' "
           f"--in_tissue_only {str(bool(_QC.get('in_tissue_only', True))).upper()}")
ALG     = config["algorithms"]

def platform_of(wc):
    """Platform of a sample: per-sample override in sample_platforms, else the global default.
    Optuna trial samples (_opt_t<N>_<sample>) resolve to their real sample name."""
    real = re.sub(r"^_opt_t\d+_", "", wc.sample)
    return config.get("sample_platforms", {}).get(real, config.get("platform", "Visium"))

# Path components never contain "/"
wildcard_constraints:
    sample = "[^/]+",
    norm   = "[^/]+",
    feat   = "[^/]+",
    alg    = "[^/]+"

ruleorder: cluster > cluster_BANKSY

import os, glob

def detect_raw_filename(wc):
    """
    Return the raw input file name (basename) in data/{sample}.
    Priority: config raw_filename if present -> .h5ad -> .h5.
    Files written by the workflow itself (01_*) are ignored, otherwise the
    chosen input would change after the first run and trigger re-runs.
    """
    d = os.path.join("data", wc.sample)
    if not os.path.isdir(d):
        return config.get("raw_filename", "filtered_feature_bc_matrix.h5")
    preferred = config.get("raw_filename", "filtered_feature_bc_matrix.h5")
    if os.path.exists(os.path.join(d, preferred)):
        return preferred
    for pat in ["*.h5ad", "*.h5"]:
        hits = sorted(h for h in glob.glob(os.path.join(d, pat))
                      if not os.path.basename(h).startswith("01_"))
        if hits:
            return os.path.basename(hits[0])
    return config.get("raw_filename", "filtered_feature_bc_matrix.h5")

# Default target: per-workflow cluster metrics for every sample listed in config["samples"]
# (norm x feature x algorithm). Summaries are built afterwards with
#   snakemake results/<sample>/clustering_metrics.csv --allowed-rules metrics_cluster_summary
rule all:
    input:
        expand("results/{sample}/{norm}/{feat}/{alg}/09_{alg}_cluster_metrics.csv",
               sample=config.get("samples", []), norm=config["normalizations"],
               feat=config["features"]["methods"], alg=list(config["algorithms"]))

# 0. Prepare input data
# ===== prepare_seurat: Seurat .rds =====
rule prepare_seurat:
    output:
        seurat = "data/{sample}/01_raw_seurat.rds"
    params:
        filename = detect_raw_filename,
        tech_platform = platform_of,
        rscript = RSCRIPT
    log: "logs/{sample}/prepare_seurat.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample} $(dirname {output.seurat}) $(dirname {log})
        echo "prepare_seurat start: {wildcards.sample} {params.tech_platform} {params.filename}" > {log} 2>&1

        # platform name is passed as written in config (Visium / VisiumHD / Stereoseq / Slideseq);
        # the scripts accept any capitalisation
        plat='{params.tech_platform}'

        fname="data/{wildcards.sample}/{params.filename}"
        if [ ! -e "$fname" ] && [ -e "{params.filename}" ]; then
            fname="{params.filename}"
        fi

        {params.rscript} scripts/prepare_seurat.R --data_dir data/{wildcards.sample} --filename "$(basename "$fname")" --out {output.seurat} --platform "$plat" {QC_ARGS} >> {log} 2>&1 || (echo "Seurat step failed, showing log:" >&2; cat {log} >&2; exit 1)
        echo "prepare_seurat done" >> {log} 2>&1
        """

# ===== prepare_anndata: anndata（raw & sc） =====
rule prepare_anndata:
    output:
        anndata = "data/{sample}/01_raw_anndata.h5ad",
        sc_anndata = "data/{sample}/01_sc_raw_anndata.h5ad"
    params:
        filename = detect_raw_filename,
        tech_platform = platform_of,
    log: "logs/{sample}/prepare_anndata.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample} $(dirname {output.anndata}) $(dirname {output.sc_anndata}) $(dirname {log})
        echo "prepare_anndata start: {wildcards.sample} {params.tech_platform} {params.filename}" > {log} 2>&1

        # platform name is passed as written in config (Visium / VisiumHD / Stereoseq / Slideseq);
        # the scripts accept any capitalisation
        plat='{params.tech_platform}'

        fname="data/{wildcards.sample}/{params.filename}"
        if [ ! -e "$fname" ] && [ -e "{params.filename}" ]; then
            fname="{params.filename}"
        fi

        # prepare_anndata.py expects filename as basename (it will join data_dir)
        fname_basename=$(basename "$fname")
        {PYTHON} scripts/prepare_anndata.py --mode raw --data_dir data/{wildcards.sample} --filename "$fname_basename" --out {output.anndata} --platform "$plat" {QC_ARGS} >> {log} 2>&1 || (echo "AnnData raw failed, showing log:" >&2; cat {log} >&2; exit 1)
        {PYTHON} scripts/prepare_anndata.py --mode sc  --data_dir data/{wildcards.sample} --filename "$fname_basename" --out {output.sc_anndata} --platform "$plat" {QC_ARGS} >> {log} 2>&1 || (echo "AnnData sc failed, showing log:" >&2; cat {log} >&2; exit 1)

        echo "prepare_anndata done" >> {log} 2>&1
        """

# 1. normalization step for LogNormalize
rule normalize_log:
    input:
        raw="data/{sample}/01_raw_seurat.rds"
    output:
        rds="data/{sample}/LogNormalize/02_LogNormalize_normalized.rds",
        plot="results/{sample}/LogNormalize/04_LogNormalize_normalize.png"
    params:
        sf=lambda wc: config["normalize"]["scale_factor"]
    log: "logs/{sample}/LogNormalize/normalize.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/LogNormalize \
                 results/{wildcards.sample}/LogNormalize \
                 logs/{wildcards.sample}/LogNormalize
        {RSCRIPT} scripts/normalize_seurat.R \
          --input {input.raw} \
          --output {output.rds} \
          --method LogNormalize \
          --scale_factor {params.sf} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# 1b. normalization step for SCTransform
rule normalize_sct:
    input:
        raw="data/{sample}/01_raw_seurat.rds"
    output:
        rds="data/{sample}/SCTransform/02_SCTransform_normalized.rds",
        hvg="results/{sample}/SCTransform/03_SCTransform_hvg.csv",
        plot="results/{sample}/SCTransform/05_SCTransform_hvg.png"
    params:
        assay=lambda wc: config["sctransform"]["assay"]
    log: "logs/{sample}/SCTransform/normalize.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/SCTransform \
                 results/{wildcards.sample}/SCTransform \
                 logs/{wildcards.sample}/SCTransform
        {RSCRIPT} scripts/sctransform_seurat.R \
          --input {input.raw} \
          --output {output.rds} \
          --assay {params.assay} \
          --hvg_csv {output.hvg} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

 # 1c. normalization step for RC (not using scale_factor)
rule normalize_rc:
     input:
         raw="data/{sample}/01_raw_seurat.rds"
     output:
         rds="data/{sample}/RC/02_RC_normalized.rds",
         plot="results/{sample}/RC/04_RC_normalize.png"
     log: "logs/{sample}/RC/normalize.log"
     shell:
         """
         mkdir -p data/{wildcards.sample}/RC \
                  results/{wildcards.sample}/RC \
                  logs/{wildcards.sample}/RC
         {RSCRIPT} scripts/normalize_seurat.R \
           --input {input.raw} \
           --output {output.rds} \
           --method RC \
           --plot_png {output.plot} \
           > {log} 2>&1
         """

 # 1d. normalization step for CLR
rule normalize_clr:
     input:
         raw="data/{sample}/01_raw_seurat.rds"
     output:
         rds="data/{sample}/CLR/02_CLR_normalized.rds",
         plot="results/{sample}/CLR/04_CLR_normalize.png"
     params:
         sf=lambda wc: config["normalize"]["scale_factor"]
     log: "logs/{sample}/CLR/normalize.log"
     shell:
         """
         mkdir -p data/{wildcards.sample}/CLR \
                  results/{wildcards.sample}/CLR \
                  logs/{wildcards.sample}/CLR
         {RSCRIPT} scripts/normalize_seurat.R \
           --input {input.raw} \
           --output {output.rds} \
           --method CLR \
           --scale_factor {params.sf} \
           --plot_png {output.plot} \
           > {log} 2>&1
         """

# 1e. normalization step for SpaNorm
rule normalize_SpaNorm:
    input:
        raw="data/{sample}/01_raw_seurat.rds"
    output:
        rds="data/{sample}/SpaNorm/02_SpaNorm_normalized.rds",
        spe="data/{sample}/SpaNorm/02_SpaNorm_spe.rds",
        plot="results/{sample}/SpaNorm/04_SpaNorm_normalize.png"
    params:
        sample_id=lambda wc: wc.sample,
        platform=platform_of
    threads: 24
    log: "logs/{sample}/SpaNorm/normalize.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/SpaNorm \
                 results/{wildcards.sample}/SpaNorm \
                 logs/{wildcards.sample}/SpaNorm
        {RSCRIPT} scripts/normalize_SpaNorm.R \
          --input {input.raw} \
          --output {output.rds} \
          --spe_output {output.spe} \
          --plot_png {output.plot} \
          --platform {params.platform} \
          --sample_id {params.sample_id} \
          --ncores {threads} \
          > {log} 2>&1
        """

# 1f. export each normalisation's values for the Python clustering methods
#     (SpaGCN / stLearn / SpatialLeiden cluster the workflow's normalised data)
rule export_normalized:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        h5="data/{sample}/{norm}/02_{norm}_normalized.h5"
    log: "logs/{sample}/{norm}/export_normalized.log"
    shell:
        """
        mkdir -p logs/{wildcards.sample}/{wildcards.norm}
        {RSCRIPT} scripts/export_normalized.R --input {input.rds} --output {output.h5} > {log} 2>&1
        """

# 0b. spot coordinates used by SpaGCN and the spatial metrics
rule make_positions:
    input:
        seurat="data/{sample}/01_raw_seurat.rds"
    output:
        csv="data/{sample}/01_positions.csv"
    log: "logs/{sample}/make_positions.log"
    shell:
        """
        mkdir -p logs/{wildcards.sample}
        {RSCRIPT} scripts/make_positions.R --input {input.seurat} --output {output.csv} > {log} 2>&1
        """

ALL_FEATURES = config["features"]["methods"]
SPECIAL = ["nnSVG", "SPARK", "SPARK-X", "mvp", "disp"]
GENERIC_FEATURES = [f for f in ALL_FEATURES if f not in SPECIAL]
FEATURE_PATTERN = "(" + "|".join(GENERIC_FEATURES) + ")"

# 2. feature selection step
rule select_features:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/{feat}/03_{feat}.rds",
        csv="results/{sample}/{norm}/{feat}/04_{feat}.csv",
        plot="results/{sample}/{norm}/{feat}/05_{feat}.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"],
        # spatial statistics (moransi / markvariogram) are scored on all detected genes unless capped
        extra=lambda wc, threads: (
            f"--n_candidates {config['features'].get('spatial_n_candidates', 0)}"
            + (f" --workers {threads}" if wc.feat == "markvariogram" else "")
        ) if wc.feat in ("moransi", "markvariogram") else ""
    threads: lambda wc: config["features"].get("markvariogram_workers", 8) if wc.feat == "markvariogram" else 1
    wildcard_constraints:
        feat=FEATURE_PATTERN
    log: "logs/{sample}/{norm}/{feat}/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat} \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat} \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}
        {RSCRIPT} scripts/select_{wildcards.feat}.R \
          --input {input.rds} \
          --output {output.rds} \
          --method {wildcards.feat} \
          --n_features {params.nfeat} \
          --csv {output.csv} \
          --plot {output.plot} \
          {params.extra} \
          > {log} 2>&1
        """

# 2c. feature selection: Mean–Variance Plot (mvp)
rule select_features_mvp:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/mvp/03_mvp.rds",
        csv="results/{sample}/{norm}/mvp/04_mvp.csv",
        plot="results/{sample}/{norm}/mvp/05_mvp.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"]
    wildcard_constraints:
        feat="mvp"
    log: "logs/{sample}/{norm}/mvp/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/mvp \
                 results/{wildcards.sample}/{wildcards.norm}/mvp \
                 logs/{wildcards.sample}/{wildcards.norm}/mvp
        {RSCRIPT} scripts/select_features.R \
          --input {input.rds} \
          --output {output.rds} \
          --method mean.var.plot \
          --n_features {params.nfeat} \
          --hvg_csv {output.csv} \
          --plot {output.plot} \
          > {log} 2>&1
        """

# 2d. feature selection: Dispersion
rule select_features_dispersion:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/disp/03_disp.rds",
        csv="results/{sample}/{norm}/disp/04_disp.csv",
        plot="results/{sample}/{norm}/disp/05_disp.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"]
    wildcard_constraints:
        feat="disp"
    log: "logs/{sample}/{norm}/disp/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/disp \
                 results/{wildcards.sample}/{wildcards.norm}/disp \
                 logs/{wildcards.sample}/{wildcards.norm}/disp
        {RSCRIPT} scripts/select_features.R \
          --input {input.rds} \
          --output {output.rds} \
          --method dispersion \
          --n_features {params.nfeat} \
          --hvg_csv {output.csv} \
          --plot {output.plot} \
          > {log} 2>&1
        """

rule select_features_nnSVG:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/nnSVG/03_nnSVG.rds",
        csv="results/{sample}/{norm}/nnSVG/04_nnSVG.csv",
        plot="results/{sample}/{norm}/nnSVG/05_nnSVG.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"],
        platform=platform_of
    wildcard_constraints:
        feat="nnSVG"
    log: "logs/{sample}/{norm}/nnSVG/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/nnSVG \
                 results/{wildcards.sample}/{wildcards.norm}/nnSVG \
                 logs/{wildcards.sample}/{wildcards.norm}/nnSVG
        {RSCRIPT} scripts/select_features_nnSVG.R \
          --input {input.rds} \
          --output {output.rds} \
          --n_features {params.nfeat} \
          --hvg_csv {output.csv} \
          --plot_png {output.plot} \
          --platform {params.platform} \
          > {log} 2>&1
        """

rule select_features_SPARK:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/SPARK/03_SPARK.rds",
        csv="results/{sample}/{norm}/SPARK/04_SPARK.csv",
        plot="results/{sample}/{norm}/SPARK/05_SPARK.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"]
    wildcard_constraints:
        feat="SPARK"
    log: "logs/{sample}/{norm}/SPARK/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/SPARK \
                 results/{wildcards.sample}/{wildcards.norm}/SPARK \
                 logs/{wildcards.sample}/{wildcards.norm}/SPARK
        {RSCRIPT} scripts/select_features_SPARK.R \
          --input {input.rds} \
          --output {output.rds} \
          --n_features {params.nfeat} \
          --hvg_csv {output.csv} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

rule select_features_SPARK_X:
    input:
        rds="data/{sample}/{norm}/02_{norm}_normalized.rds"
    output:
        rds="data/{sample}/{norm}/SPARK-X/03_SPARK-X.rds",
        csv="results/{sample}/{norm}/SPARK-X/04_SPARK-X.csv",
        plot="results/{sample}/{norm}/SPARK-X/05_SPARK-X.png"
    params:
        nfeat=lambda wc: config["features"]["n_features"]
    wildcard_constraints:
        feat="SPARK-X"
    log: "logs/{sample}/{norm}/SPARK-X/select.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/SPARK-X \
                 results/{wildcards.sample}/{wildcards.norm}/SPARK-X \
                 logs/{wildcards.sample}/{wildcards.norm}/SPARK-X
        {RSCRIPT} scripts/select_features_SPARK-X.R \
          --input {input.rds} \
          --output {output.rds} \
          --n_features {params.nfeat} \
          --hvg_csv {output.csv} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# 3. dimensionality reduction step
rule reduce:
    input:
        rds="data/{sample}/{norm}/{feat}/03_{feat}.rds"
    output:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds",
        embeddings="results/{sample}/{norm}/{feat}/05_umap.csv",
        plot="results/{sample}/{norm}/{feat}/06_umap.png"
    params:
        npcs=lambda wc: step3_npcs(),
        umap_dims=lambda wc: n_pcs("clustering"),
        platform=platform_of
    log: "logs/{sample}/{norm}/{feat}/reduce.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat} \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat} \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}
        {RSCRIPT} scripts/reduce_dimensions.R \
          --input {input.rds} \
          --output {output.rds} \
          --npcs {params.npcs} \
          --umap_dims {params.umap_dims} \
          --umap_csv {output.embeddings} \
          --plot_png {output.plot} \
          --platform {params.platform} \
          > {log} 2>&1
        """

# for clustering algorithms outside louvain/leiden within 'rule cluster', they will run corresponding script such as rule 'cluster_Banksy'
GENERIC = [a for a in ALG.keys()
           if a not in ("BANKSY", "BayesSpace", "SpaGCN", "stlearn", "SpatialPCA", "SpatialLeiden")]
GENERIC_PATTERN = "(" + "|".join(GENERIC) + ")"

# 4.1 cluster inside Seurat
rule cluster:
    input:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds"
    output:
        rds ="data/{sample}/{norm}/{feat}/{alg}/06_{alg}_clustered.rds",
        csv ="results/{sample}/{norm}/{feat}/{alg}/07_{alg}_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/{alg}/08_{alg}_clusters.png"
    params:
        res=lambda wc: config["clustering"]["resolution"],
        algnum=lambda wc: ALG[wc.alg],
        dims=lambda wc: n_pcs("clustering")
    wildcard_constraints:
        alg=GENERIC_PATTERN
    log: "logs/{sample}/{norm}/{feat}/{alg}/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg} \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg} \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg}
        CLM=""
        if [ "{wildcards.alg}" = "leiden" ]; then
            CLM="--clustermethod igraph"
        fi
        {RSCRIPT} scripts/cluster_seurat.R \
          --input  {input.rds} \
          --output {output.rds} \
          --resolution {params.res} \
          --algorithm {params.algnum} \
          --dims {params.dims} \
          $CLM \
          --meta_csv {output.csv} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# 4.2. BANKSY‐specific clustering rule
rule cluster_BANKSY:
    input:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds"
    output:
        rds ="data/{sample}/{norm}/{feat}/BANKSY/06_BANKSY_clustered.rds",
        csv ="results/{sample}/{norm}/{feat}/BANKSY/07_BANKSY_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/BANKSY/08_BANKSY_clusters.png"
    params:
        res = lambda wc: config["BANKSY"]["resolution"],
        lam = lambda wc: config["BANKSY"]["lambda"],
        k_g = lambda wc: config["BANKSY"]["k_geom"],
        extra = lambda wc: " ".join(filter(None, [
            f"--spatial_mode {config['BANKSY'].get('spatial_mode', 'kNN_median')}",
            f"--npcs {n_pcs('BANKSY')}",
            f"--k_neighbors {config['BANKSY'].get('k_neighbors', 50)}",
            "--use_agf" if config['BANKSY'].get('use_agf', False) else "",
            f"--M {config['BANKSY']['M']}" if config['BANKSY'].get('M') is not None else ""]))
    log: "logs/{sample}/{norm}/{feat}/BANKSY/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BANKSY \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BANKSY \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BANKSY

        {RSCRIPT} scripts/cluster_Banksy.R \
          --input  {input.rds} \
          --output {output.rds} \
          --resolution {params.res} \
          --lambda {params.lam} \
          --k_geom {params.k_g} \
          {params.extra} \
          --meta_csv {output.csv} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# 4.3. BayesSpace‐specific clustering rule
rule cluster_BayesSpace:
    input:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds"
    output:
        csv ="results/{sample}/{norm}/{feat}/BayesSpace/07_BayesSpace_clusters.csv",
        rds="data/{sample}/{norm}/{feat}/BayesSpace/06_BayesSpace_clustered.rds"
    params:
        sce_rds="data/{sample}/{norm}/{feat}/BayesSpace/05_BayesSpace_sce.rds",
        meta_csv="results/{sample}/{norm}/{feat}/BayesSpace/07_BayesSpace_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/BayesSpace/08_BayesSpace_clusters.png",
        platform=platform_of,
        d=lambda wc: n_pcs("BayesSpace")
    log: "logs/{sample}/{norm}/{feat}/BayesSpace/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BayesSpace \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BayesSpace \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/BayesSpace

        {RSCRIPT} scripts/cluster_BayesSpace.R \
          --input {input.rds} \
          --output {output.rds} \
          --sce_output {params.sce_rds} \
          --meta_csv {params.meta_csv} \
          --plot_png {params.plot} \
          --q {config[BayesSpace][q]} \
          --gamma {config[BayesSpace][gamma]} \
          --nrep {config[BayesSpace][nrep]} \
          --burn_in {config[BayesSpace][burn_in]} \
          --d {params.d} \
          --model {config[BayesSpace][model]} \
          --init_method {config[BayesSpace][init_method]} \
          --precision {config[BayesSpace][precision]} \
          --platform {params.platform} \
          > {log} 2>&1
        """

# 4.4 SpaGCN‐specific clustering rule
rule cluster_SpaGCN:
    input:
        raw_h5ad     = "data/{sample}/01_raw_anndata.h5ad",
        positions    = "data/{sample}/01_positions.csv",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv",
        norm_h5      = "data/{sample}/{norm}/02_{norm}_normalized.h5",
        # Ensure image is an input if you use histology!
        img = lambda wc: f"data/{wc.sample}/histology.tif" if config["SpaGCN"].get("histology", False) else []
    output:
        rds          = "data/{sample}/{norm}/{feat}/SpaGCN/06_SpaGCN_clustered.rds",
        csv          = "results/{sample}/{norm}/{feat}/SpaGCN/07_SpaGCN_clusters.csv",
        png          = "results/{sample}/{norm}/{feat}/SpaGCN/08_SpaGCN_clusters.png"
    params:
        refine_flag  = lambda wc: "--refine" if config["SpaGCN"]["refine"] else "",
        hist_flag    = lambda wc: "--histology" if config["SpaGCN"].get("histology", False) else "",
        img_cmd      = lambda wc, input: f"--img {input.img[0]}" if input.img else "",
        beta         = config["SpaGCN"]["beta"],
        alpha        = config["SpaGCN"]["alpha"],
        n_clusters   = config["SpaGCN"]["n_clusters"],
        num_pcs      = lambda wc: n_pcs("SpaGCN"),
        platform     = platform_of,
        env          = ENV_SPAGCN
    threads: 36
    log: "logs/{sample}/{norm}/{feat}/SpaGCN/cluster.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN
        
        # requires micromamba on PATH (e.g. `module load micromamba` on HPC)
        micromamba run -p {params.env} python scripts/cluster_SpaGCN.py \
          --raw_h5ad   {input.raw_h5ad} \
          --positions  {input.positions} \
          --features   {input.features_csv} \
          --norm_h5    {input.norm_h5} \
          --out_rds    {output.rds} \
          --out_csv    {output.csv} \
          --out_png    {output.png} \
          {params.img_cmd} \
          --n_clusters {params.n_clusters} \
          --num_pcs    {params.num_pcs} \
          --platform   {params.platform} \
          --beta       {params.beta} \
          --alpha      {params.alpha} \
          {params.refine_flag} {params.hist_flag} \
          > {log} 2>&1
        """

# 4.5. stLearn–specific clustering rule
rule cluster_stlearn:
    input:
        raw_h5ad     = "data/{sample}/01_sc_raw_anndata.h5ad",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv",
        norm_h5      = "data/{sample}/{norm}/02_{norm}_normalized.h5"
    output:
        h5ad        = "data/{sample}/{norm}/{feat}/stlearn/09_stlearn_clustered.h5ad",
        csv         = "results/{sample}/{norm}/{feat}/stlearn/07_stlearn_clusters.csv",
        png         = "results/{sample}/{norm}/{feat}/stlearn/08_stlearn_clusters.png",
        rds         = "data/{sample}/{norm}/{feat}/stlearn/06_stlearn_clustered.rds"
    params:
        morph_flag  = "--use-morphology" if config["stlearn"]["use_morphology"] else "",
        n_pcs       = lambda wc: n_pcs("stlearn"),
        resolution  = lambda wc: config["stlearn"].get("resolution", 1.0),
        n_neighbors = lambda wc: config["stlearn"]["n_neighbors"],
        platform    = platform_of,
        env         = ENV_STLEARN,
        extra       = lambda wc: " ".join(filter(None, [
            f"--sme-weights {config['stlearn'].get('sme_weights', 'weights_matrix_all')}",
            f"--clustering-method {config['stlearn'].get('clustering_method', 'louvain')}",
            f"--n-clusters {config['stlearn']['n_clusters']}" if config['stlearn'].get('n_clusters') is not None else "",
            "--use-localization" if config['stlearn'].get('use_localization', False) else "",
            f"--localization-eps {config['stlearn'].get('localization_eps', 20.0)}",
            f"--localization-min-samples {config['stlearn'].get('localization_min_samples', 1)}"]))
    threads: 36
    log: "logs/{sample}/{norm}/{feat}/stlearn/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn
        
        # requires micromamba on PATH (e.g. `module load micromamba` on HPC)
        micromamba run -p {params.env} python scripts/cluster_stlearn.py \
          --raw-h5ad   {input.raw_h5ad} \
          --features   {input.features_csv} \
          --norm-h5    {input.norm_h5} \
          --out-h5ad   {output.h5ad} \
          --out-csv    {output.csv} \
          --out-rds    {output.rds} \
          --out-png    {output.png} \
          --platform   {params.platform} \
          {params.morph_flag} \
          --n-pcs      {params.n_pcs} \
          --resolution {params.resolution} \
          --n-neighbors {params.n_neighbors} \
          {params.extra} \
          > {log} 2>&1
        """

# 4.6 SpatialPCA‐specific clustering rule
rule cluster_SpatialPCA:
    input:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds"
    output:
        rds="data/{sample}/{norm}/{feat}/SpatialPCA/06_SpatialPCA_clustered.rds",
        csv="results/{sample}/{norm}/{feat}/SpatialPCA/07_SpatialPCA_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/SpatialPCA/08_SpatialPCA_clusters.png"
    params:
        n_clusters=lambda wc: config["SpatialPCA"]["n_clusters"],
        gene_number=lambda wc: config["features"]["n_features"],
        refine_flag=lambda wc: f"--refine {str(config['SpatialPCA']['refine']).upper()}",
        spatial_pcs=lambda wc: n_pcs("SpatialPCA")
    log: "logs/{sample}/{norm}/{feat}/SpatialPCA/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialPCA \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialPCA \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialPCA

        {RSCRIPT} scripts/cluster_SpatialPCA.R \
          --input      {input.rds} \
          --output     {output.rds} \
          --meta_csv   {output.csv} \
          --plot_png   {output.plot} \
          --n_clusters {params.n_clusters} \
          --gene_number {params.gene_number} \
          --spatial_pcs {params.spatial_pcs} \
          {params.refine_flag} \
          > {log} 2>&1
        """

# 4.6 SpatialLeiden‐specific clustering rule
rule cluster_SpatialLeiden:
    input:
        raw_h5ad     = "data/{sample}/01_sc_raw_anndata.h5ad",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv",
        norm_h5      = "data/{sample}/{norm}/02_{norm}_normalized.h5"
    output:
        h5ad = "data/{sample}/{norm}/{feat}/SpatialLeiden/09_SpatialLeiden_clustered.h5ad",
        csv  = "results/{sample}/{norm}/{feat}/SpatialLeiden/07_SpatialLeiden_clusters.csv",
        png  = "results/{sample}/{norm}/{feat}/SpatialLeiden/08_SpatialLeiden_clusters.png",
        rds  = "data/{sample}/{norm}/{feat}/SpatialLeiden/06_SpatialLeiden_clustered.rds"
    params:
        n_neighbors = lambda wc: config["SpatialLeiden"]["n_neighbors"],
        n_pcs      = lambda wc: n_pcs("SpatialLeiden"),
        resolution  = lambda wc: config["SpatialLeiden"]["resolution"],
        layer_ratio = lambda wc: config["SpatialLeiden"]["layer_ratio"],
        coord_type  = lambda wc: config["SpatialLeiden"]["coord_type"],
        refine_flag = lambda wc: "--refine" if config["SpatialLeiden"]["refine"] else "",
        python      = PY_SPATIALLEIDEN,
        extra       = lambda wc: " ".join([
            "--use-weights " + " ".join(str(v).lower() for v in config["SpatialLeiden"].get("use_weights", [True, True])),
            "--search-start " + " ".join(str(v) for v in config["SpatialLeiden"].get("search_start", [1.0, 0.4])),
            f"--n-leiden-iter {config['SpatialLeiden'].get('n_leiden_iter', -1)}"]),
        search_n_clusters_flag = lambda wc: (
            f"--search-n-clusters {config['SpatialLeiden']['search_n_clusters']}"
        ) if config["SpatialLeiden"].get("search_n_clusters", None) is not None else ""
    threads: 36
    log: "logs/{sample}/{norm}/{feat}/SpatialLeiden/cluster.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialLeiden \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialLeiden \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpatialLeiden
        
        {params.python} scripts/cluster_SpatialLeiden.py \
          --raw-h5ad {input.raw_h5ad} \
          --features {input.features_csv} \
          --norm-h5 {input.norm_h5} \
          --out-h5ad {output.h5ad} \
          --out-csv {output.csv} \
          --out-png {output.png} \
          --out-rds {output.rds} \
          --n-neighbors {params.n_neighbors} \
          --n-pcs {params.n_pcs} \
          --resolution {params.resolution} \
          --layer-ratio {params.layer_ratio} \
          --coord-type {params.coord_type} \
          {params.search_n_clusters_flag} \
          {params.refine_flag} \
          {params.extra} \
          > {log} 2>&1
        """

# 5.1. Cell type annotation & metrics computation
import glob, os
# --- annotation step ---
rule annotate_clusters:
    input:
        seurat   = "data/{sample}/{norm}/{feat}/{alg}/06_{alg}_clustered.rds",
        clusters = "results/{sample}/{norm}/{feat}/{alg}/07_{alg}_clusters.csv",
        truth    = "data/Ground_truth/{sample}/labels.tsv"
    output:
        seurat_annot = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_annotated.rds",
        seurat_plot  = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_annotation.png"
    log: "logs/{sample}/{norm}/{feat}/{alg}/annotate.log"
    params:
        save_rds = config.get("save_annotation_rds", True)
    shell:
        """
        mkdir -p results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg}
        {RSCRIPT} scripts/annotate_clusters.R \
          --seurat      {input.seurat} \
          --clusters    {input.clusters} \
          --truth       {input.truth} \
          --seurat_out  {output.seurat_annot} \
          --plot_seurat {output.seurat_plot} \
          --save_rds    {params.save_rds} \
          > {log} 2>&1
        """

# --- metrics computation step ---
rule compute_metrics:
    input:
        seurat_annot = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_annotated.rds",
        truth        = "data/Ground_truth/{sample}/labels.tsv"
    output:
        metrics = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_metrics.csv",
        plot    = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_metrics_plot.png"
    log:
        "logs/{sample}/{norm}/{feat}/{alg}/metrics.log"
    shell:
        """
        mkdir -p results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg}
        {RSCRIPT} scripts/compute_metrics.R \
          --seurat_annot {input.seurat_annot} \
          --truth         {input.truth} \
          --out           {output.metrics} \
          --plot          {output.plot} \
          --algorithm     {wildcards.alg} \
          > {log} 2>&1
        """

# 5.2. summary all metrics & simple plotting
rule metrics_summary:
    input:
        # Dynamically collect the already generated metrics files, with recursive=True to search through subdirectories
        metrics=lambda wc: [f for f in glob.glob(
            os.path.join("results", wc.sample, "**", "09_*_metrics.csv"),
            recursive=True) if not f.endswith("_cluster_metrics.csv")]
    output:
        summary="results/{sample}/metrics/metrics_summary.csv",
        plot   ="results/{sample}/metrics/metrics_summary.png"
    params:
        sample=lambda wc: wc.sample
    log:
        "logs/{sample}/metrics_summary.log"
    run:
        # if not found, stop
        if not input.metrics:
            raise ValueError(f"No metrics CSV found for sample {params.sample}")
        shell("mkdir -p results/{params.sample}/metrics")
        mf = " ".join(input.metrics)
        shell("mkdir -p results/{params.sample}/metrics")
        shell(
            f"{RSCRIPT} scripts/aggregate_metrics.R "
            f"--metrics_dir results/{params.sample} "
            f"--out_csv {output.summary} "
            f"--out_plot {output.plot} "
            f"> {log} 2>&1"
        )

# 5.3. Cluster-only metrics (no annotation step; raw cluster IDs vs truth)
# These avoid the data-leakage bias introduced when SingleR is trained on
# the same ground-truth labels used for evaluation.
# Samples whose ground-truth file only annotates a subset of spots;
# everything not listed is filled with "_Other_" so NMI/ARI/HOM/COM are
# well-defined on the full clustered set.
PARTIAL_TRUTH_SAMPLES = {
    "Mouse_embryo_Visium",
    "Mouse_embryo_E12_VisiumHD",
}

rule compute_cluster_metrics:
    input:
        clusters  = "results/{sample}/{norm}/{feat}/{alg}/07_{alg}_clusters.csv",
        positions = "data/{sample}/01_positions.csv",
        truth     = "data/Ground_truth/{sample}/labels.tsv"
    output:
        metrics = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_cluster_metrics.csv",
        plot    = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_cluster_metrics_plot.png"
    log:
        "logs/{sample}/{norm}/{feat}/{alg}/cluster_metrics.log"
    params:
        gt_label = lambda wc: "partial" if re.sub(r"^_opt_t\d+_", "", wc.sample) in PARTIAL_TRUTH_SAMPLES else "full"
    shell:
        """
        mkdir -p results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg}
        {RSCRIPT} scripts/compute_cluster_metrics.R \
          --clusters            {input.clusters} \
          --positions           {input.positions} \
          --truth               {input.truth} \
          --out                 {output.metrics} \
          --plot                {output.plot} \
          --algorithm           {wildcards.alg} \
          --ground_truth_label  {params.gt_label} \
          > {log} 2>&1
        """

rule metrics_cluster_summary:
    input:
        metrics = lambda wc: glob.glob(
            os.path.join("results", wc.sample, "**", "09_*_cluster_metrics.csv"),
            recursive=True
        )
    output:
        summary = "results/{sample}/clustering_metrics.csv"
    params:
        sample = lambda wc: wc.sample
    log:
        "logs/{sample}/metrics_cluster_summary.log"
    run:
        if not input.metrics:
            raise ValueError(f"No cluster metrics CSV found for sample {params.sample}")
        shell("mkdir -p results/{params.sample}")
        shell(
            f"{RSCRIPT} scripts/aggregate_cluster_metrics.R "
            f"--metrics_dir results/{params.sample} "
            f"--out_csv {output.summary} "
            f"> {log} 2>&1"
        )
