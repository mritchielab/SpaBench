# SpaBench
<img  src="logo/Snakemake.png">
Using Snakemake, SpaBench is a comprehensive spatial transcriptomics benchmarking workflow which implements 23 distinct methods (contain non-spatial and spatially-aware ones) and includes 360 different combinations of workflows. SpaBench begins from raw spatially-resolved gene count matrix, followed by normalization, feature selection, dimensionality reduction, and clustering. As the end point, SpaBench utilizes SingleR for cell type / spatial domain annotation and computes QC metrics, including both non-spatial ones such as NMI, HOM, COM, CHAOS, PAS, ASW, and spatially-aware IISC (Integrated Image–Spatial Cohesion).

## To run the whole workflow:
```
module load snakemake
snakemake -j 16
```

## To run specific branch of workflow:
```
snakemake -j 16 -p   /pathtodata/{CLR,LogNormalize,RC,SCTransform}/{vst,mvp,disp,markvariogram,moransi,nnSVG,SPARK,SPARK-X}/SpaGCN/06_SpaGCN_clustered.rds
```

## To run particular steps:
```
snakemake data/CLR/nnSVG/BayesSpace/06_BayesSpace_clustered.rds -j 16 -p --allowed-rules cluster_BayesSpace 
```

## To run cell type / spatial domain annotation with QC metrics computation:
```
snakemake   $(find results/data -type f -name "07_*_clusters.csv" \
     | sed -e 's#/07_\(.*\)_clusters.csv#/09_\1_metrics.csv#')   --cores 16 -p --allowed-rules compute_metrics
```

## Resource Management
Snakemake can automatically determine which parts of the workflow can be run in parallel. By specifying more than one available core:
```
snakemake --cores 16 /pathtodata/{CLR,LogNormalize,RC,SCTransform}/{vst,mvp,disp,markvariogram,moransi,nnSVG,SPARK,SPARK-X}/{BANKSY, BayesSpace. leiden, louvain, louvain_refine, SLM, SpatialPCA, stlearn, SpaGCN}/06_*_clustered.rds
```
