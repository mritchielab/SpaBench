# SpaBench
<img  src="logo/Snakemake.png">
Building upon Snakemake, SpaBench is a comprehensive spatial transcriptomics benchmarking workflow which implements 23 distinct methods (contain non-spatial and spatially-aware ones) and includes 360 different combinations of workflows. SpaBench begins from raw spatially-resolved gene count matrix, followed by normalization, feature selection, dimensionality reduction, and clustering. As the end point, SpaBench utilizes SingleR for cell type / spatial domain annotation and computes QC metrics, including both non-spatial ones such as NMI, HOM, COM, CHAOS, PAS, ASW, and spatially-aware IISC (Integrated Image–Spatial Cohesion).

## Overview of Algorithm Implementations
| Algorithm | Functional Step | Language | Citation |
| --- | --- | --- | --- |
| `Seurat` | Object construction & Normalization & Feature selection & Dimensionality reduction & Clustering | R | [Hao et al. (2023)](https://www.nature.com/articles/s41587-023-01767-y)
| `Squidpy` | Object construction & Image feature computation | Python | [Palla et al. (2022)](https://www.nature.com/articles/s41592-021-01358-2)
| `SpatialExperiment` | Object construction | R | [Righelli et al. (2022)](https://academic.oup.com/bioinformatics/article/38/11/3128/6575443)
| `SpaNorm` | Normalization | R | [Salim et al. (2025)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03565-y)
| `SPARK` & `SPARK-X` | Feature selection | R | [Weber et al. (2023)](https://www.nature.com/articles/s41467-023-39748-z)
| `nnSVG` | Feature selection | R | [Zhu et al. (2021)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02404-0)
| `SpatialPCA` |  Dimensionality reduction & Clustering | R | [Shang et al. (2022)](https://www.nature.com/articles/s41467-022-34879-1)
| `BANKSY` | Clustering | R & Python | [Singhal et al. (2023)](https://www.nature.com/articles/s41588-024-01664-3)
| `stlearn` | Clustering | Python | [Pham et al. (2023)](https://www.nature.com/articles/s41467-023-43120-6 )
| `BayesSpace` | Clustering | R | [Zhao et al. (2021)](https://www.nature.com/articles/s41587-021-00935-2)
| `SpaGCN` | Clustering | Python | [Hu et al. (2021)](https://www.nature.com/articles/s41592-021-01255-8)

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
## Before running this workflow:
### Input data format
Named as ```01_raw_seurat.rds``` under ```data/samplename/```
