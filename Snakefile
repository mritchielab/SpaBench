# Snakefile

configfile: "config.yaml"

RSCRIPT = "/stornext/System/data/software/rhel/9/base/tools/R/4.5.1/lib64/R/bin/Rscript"
SAMPLES = config["samples"]
ALG     = config["algorithms"]

ruleorder: cluster > cluster_BANKSY

# 0. Prepare input data
# ===== prepare_seurat: Seurat .rds =====
rule prepare_seurat:
    output:
        seurat = "data/{sample}/01_raw_seurat.rds"
    params:
        filename = detect_raw_filename,
        tech_platform = lambda wc: config.get("sample_platforms", {}).get(wc.sample, config.get("platform", "Visium")),
        rscript = RSCRIPT
    log: "logs/{sample}/prepare_seurat.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample} $(dirname {output.seurat}) $(dirname {log})
        echo "prepare_seurat start: {wildcards.sample} {params.tech_platform} {params.filename}" > {log} 2>&1

        plat=$(echo '{params.tech_platform}' | tr '[:upper:]' '[:lower:]')
        case "$plat" in
            stereoseq) plat="StereoSeq" ;;
            slideseq) plat="SlideSeq" ;;
            visium*|visiumhd) plat="Visium" ;;
            *) plat="{params.tech_platform}" ;;
        esac

        fname="data/{wildcards.sample}/{params.filename}"
        if [ ! -e "$fname" ] && [ -e "{params.filename}" ]; then
            fname="{params.filename}"
        fi

        {params.rscript} scripts/prepare_seurat.R --data_dir data/{wildcards.sample} --filename "$fname" --out {output.seurat} --platform "$plat" >> {log} 2>&1 || (echo "Seurat step failed, showing log:" >&2; cat {log} >&2; exit 1)
        echo "prepare_seurat done" >> {log} 2>&1
        """

# ===== prepare_anndata: anndata（raw & sc） =====
rule prepare_anndata:
    output:
        anndata = "data/{sample}/01_raw_anndata.h5ad",
        sc_anndata = "data/{sample}/01_sc_raw_anndata.h5ad"
    params:
        filename = detect_raw_filename,
        tech_platform = lambda wc: config.get("sample_platforms", {}).get(wc.sample, config.get("platform", "Visium")),
    log: "logs/{sample}/prepare_anndata.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample} $(dirname {output.anndata}) $(dirname {output.sc_anndata}) $(dirname {log})
        echo "prepare_anndata start: {wildcards.sample} {params.tech_platform} {params.filename}" > {log} 2>&1

        plat=$(echo '{params.tech_platform}' | tr '[:upper:]' '[:lower:]')
        case "$plat" in
            stereoseq) plat="StereoSeq" ;;
            slideseq) plat="SlideSeq" ;;
            visium*|visiumhd) plat="Visium" ;;
            *) plat="{params.tech_platform}" ;;
        esac

        fname="data/{wildcards.sample}/{params.filename}"
        if [ ! -e "$fname" ] && [ -e "{params.filename}" ]; then
            fname="{params.filename}"
        fi

        # prepare_anndata.py expects filename as basename (it will join data_dir)
        fname_basename=$(basename "$fname")
        python scripts/prepare_anndata.py --mode raw --data_dir data/{wildcards.sample} --filename "$fname_basename" --out {output.anndata} --platform "$plat" >> {log} 2>&1 || (echo "AnnData raw failed, showing log:" >&2; cat {log} >&2; exit 1)
        python scripts/prepare_anndata.py --mode sc  --data_dir data/{wildcards.sample} --filename "$fname_basename" --out {output.sc_anndata} --platform "$plat" >> {log} 2>&1 || (echo "AnnData sc failed, showing log:" >&2; cat {log} >&2; exit 1)

        echo "prepare_anndata done" >> {log} 2>&1
        """

# 1. normalization step for LogNormalize
rule normalize_log:
    input:
        raw="data/{sample}/01_raw_seurat.rds"
    output:
        rds="data/{sample}/LogNormalize/02_LogNormalize_normalized.rds",
        counts="results/{sample}/LogNormalize/03_LogNormalize_counts.csv",
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
          --counts_csv {output.counts} \
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
        scaled="results/{sample}/SCTransform/04_SCTransform_scaled.csv",
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
          --scaled_csv {output.scaled} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

 # 1c. normalization step for RC (not using scale_factor)
rule normalize_rc:
     input:
         raw="data/{sample}/01_raw_seurat.rds"
     output:
         rds="data/{sample}/RC/02_RC_normalized.rds",
         counts="results/{sample}/RC/03_RC_counts.csv",
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
           --counts_csv {output.counts} \
           --plot_png {output.plot} \
           > {log} 2>&1
         """

 # 1d. normalization step for CLR
rule normalize_clr:
     input:
         raw="data/{sample}/01_raw_seurat.rds"
     output:
         rds="data/{sample}/CLR/02_CLR_normalized.rds",
         counts="results/{sample}/CLR/03_CLR_counts.csv",
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
           --counts_csv {output.counts} \
           --plot_png {output.plot} \
           > {log} 2>&1
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
        nfeat=lambda wc: config["features"]["n_features"]
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
        platform=lambda wc: config["platform"]
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
        npcs=lambda wc: config["pca"]["npcs"]
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
          --umap_csv {output.embeddings} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# for clustering algorithms outside louvain/leiden within 'rule cluster', they will run corresponding script such as rule 'cluster_Banksy'
GENERIC = [a for a in ALG.keys()
           if a not in ("BANKSY", "BayesSpace", "SpaGCN", "stlearn", "SpatialPCA")]
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
        algnum=lambda wc: ALG[wc.alg]
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
        res=lambda wc: config["clustering"]["resolution"]
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
          --meta_csv {output.csv} \
          --plot_png {output.plot} \
          > {log} 2>&1
        """

# 4.3. BayesSpace‐specific clustering rule
rule cluster_BayesSpace:
    input:
        rds="data/{sample}/{norm}/{feat}/04_reduced.rds"
    output:
        rds="data/{sample}/{norm}/{feat}/BayesSpace/06_BayesSpace_clustered.rds",
        csv ="results/{sample}/{norm}/{feat}/BayesSpace/07_BayesSpace_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/BayesSpace/08_BayesSpace_clusters.png"
    params:
        sce_rds="data/{sample}/{norm}/{feat}/BayesSpace/05_BayesSpace_sce.rds",
        meta_csv="results/{sample}/{norm}/{feat}/BayesSpace/07_BayesSpace_clusters.csv",
        plot="results/{sample}/{norm}/{feat}/BayesSpace/08_BayesSpace_clusters.png"
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
          --q {config[clustering][q]} \
          --nrep {config[clustering][nrep]} \
          --burn_in {config[clustering][burn_in]} \
          > {log} 2>&1
        """

# 4.4 SpaGCN‐specific clustering rule
rule cluster_SpaGCN:
    input:
        raw_h5ad     = "data/{sample}/01_raw_anndata.h5ad",
        positions    = "data/{sample}/01_positions.csv",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv"
    output:
        rds         = "data/{sample}/{norm}/{feat}/SpaGCN/06_SpaGCN_clustered.rds",
        csv         = "results/{sample}/{norm}/{feat}/SpaGCN/07_SpaGCN_clusters.csv",
        png         = "results/{sample}/{norm}/{feat}/SpaGCN/08_SpaGCN_clusters.png"
    params:
        refine_flag=lambda wc: "--refine" if config["SpaGCN"]["refine"] else ""
    threads: 36
    log: "logs/{sample}/{norm}/{feat}/SpaGCN/cluster.log"
    shell:
        r"""
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/SpaGCN
        
        eval "$(micromamba shell hook --shell bash)"
        module load micromamba

        micromamba run -n SpaGCN-env python scripts/cluster_SpaGCN.py \
          --raw_h5ad   {input.raw_h5ad} \
          --positions  {input.positions} \
          --features   {input.features_csv} \
          --out_rds    {output.rds} \
          --out_csv    {output.csv} \
          --out_png    {output.png} \
          {params.refine_flag} \
          > {log} 2>&1
        """

# 4.5. stLearn–specific clustering rule
rule cluster_stlearn:
    input:
        raw_h5ad     = "data/{sample}/01_sc_raw_anndata.h5ad",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv"
    output:
        h5ad        = "data/{sample}/{norm}/{feat}/stlearn/09_stlearn_clustered.h5ad",
        csv         = "results/{sample}/{norm}/{feat}/stlearn/07_stlearn_clusters.csv",
        png         = "results/{sample}/{norm}/{feat}/stlearn/08_stlearn_clusters.png",
        rds         = "data/{sample}/{norm}/{feat}/stlearn/06_stlearn_clustered.rds"
    params:
        morph_flag  = "--use-morphology" if config["stlearn"]["use_morphology"] else "",
        n_pcs       = lambda wc: config["stlearn"].get("n_pcs", 50),
        resolution  = lambda wc: config["stlearn"].get("resolution", 1.0),
        n_neighbors = lambda wc: config["stlearn"]["n_neighbors"],
        platform    = lambda wc: config["platform"]
    threads: 36
    log: "logs/{sample}/{norm}/{feat}/stlearn/cluster.log"
    shell:
        """
        mkdir -p data/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn \
                 results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn \
                 logs/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/stlearn
        
        eval "$(micromamba shell hook --shell bash)"
        module load micromamba

        micromamba run -n stlearn-env python scripts/cluster_stlearn.py \
          --raw-h5ad   {input.raw_h5ad} \
          --features   {input.features_csv} \
          --out-h5ad   {output.h5ad} \
          --out-csv    {output.csv} \
          --out-rds    {output.rds} \
          --out-png    {output.png} \
          --platform   {params.platform} \
          {params.morph_flag} \
          --n-pcs      {params.n_pcs} \
          --resolution {params.resolution} \
          --n-neighbors {params.n_neighbors} \
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
        refine_flag=lambda wc: f"--refine {str(config['SpatialPCA']['refine']).upper()}"
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
          {params.refine_flag} \
          > {log} 2>&1
        """

# 4.7 SpatialLeiden‐specific clustering rule
rule cluster_SpatialLeiden:
    input:
        raw_h5ad     = "data/{sample}/01_sc_raw_anndata.h5ad",
        features_csv = "results/{sample}/{norm}/{feat}/04_{feat}.csv"
    output:
        h5ad = "data/{sample}/{norm}/{feat}/SpatialLeiden/09_SpatialLeiden_clustered.h5ad",
        csv  = "results/{sample}/{norm}/{feat}/SpatialLeiden/07_SpatialLeiden_clusters.csv",
        png  = "results/{sample}/{norm}/{feat}/SpatialLeiden/08_SpatialLeiden_clusters.png",
        rds  = "data/{sample}/{norm}/{feat}/SpatialLeiden/06_SpatialLeiden_clustered.rds"
    params:
        n_neighbors = lambda wc: config["SpatialLeiden"]["n_neighbors"],
        n_pcs      = lambda wc: config["SpatialLeiden"]["n_pcs"],
        resolution  = lambda wc: config["SpatialLeiden"]["resolution"],
        layer_ratio = lambda wc: config["SpatialLeiden"]["layer_ratio"],
        coord_type  = lambda wc: config["SpatialLeiden"]["coord_type"],
        spatial_n_neighbors = lambda wc: config["SpatialLeiden"]["spatial_n_neighbors"],
        refine_flag = lambda wc: "--refine" if config["SpatialLeiden"]["refine"] else "",
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
        
        module load python
        python scripts/cluster_SpatialLeiden.py \
          --raw-h5ad {input.raw_h5ad} \
          --features {input.features_csv} \
          --out-h5ad {output.h5ad} \
          --out-csv {output.csv} \
          --out-png {output.png} \
          --out-rds {output.rds} \
          --n-neighbors {params.n_neighbors} \
          --n-pcs {params.n_pcs} \
          --resolution {params.resolution} \
          --layer-ratio {params.layer_ratio} \
          --coord-type {params.coord_type} \
          --spatial-n-neighs {params.spatial_n_neighbors} \
          {params.search_n_clusters_flag} \
          {params.refine_flag} \
          > {log} 2>&1
        """

# 5.1. Cell type annotation & metrics computation
import glob, os
rule compute_metrics:
    input:
        seurat   = "data/{sample}/{norm}/{feat}/{alg}/06_{alg}_clustered.rds",
        clusters = "results/{sample}/{norm}/{feat}/{alg}/07_{alg}_clusters.csv",
        truth    = "data/Ground_truth/{sample}/labels.tsv"
    output:
        metrics  = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_metrics.csv",
        seurat_annot = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_annotated.rds",
        seurat_plot  = "results/{sample}/{norm}/{feat}/{alg}/09_{alg}_annotation.png"
    #threads: 1
    #resources:
      #mem_mb=16000 
    log:
        "logs/{sample}/{norm}/{feat}/{alg}/metrics.log"
    shell:
        """
        mkdir -p results/{wildcards.sample}/{wildcards.norm}/{wildcards.feat}/{wildcards.alg}
        {RSCRIPT} scripts/compute_metrics.R \
          --seurat   {input.seurat} \
          --clusters {input.clusters} \
          --truth    {input.truth} \
          --out      {output.metrics} \
          --seurat_out   {output.seurat_annot} \
          --plot_seurat  {output.seurat_plot} \
          > {log} 2>&1
        """

# 5.2. summary all metrics & simple plotting
rule metrics_summary:
    input:
        # Dynamically collect the already generated metrics files, with recursive=True to search through subdirectories
        metrics=lambda wc: glob.glob(
            os.path.join("results", wc.sample, "**", "09_*_metrics.csv"),
            recursive=True
        )
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
