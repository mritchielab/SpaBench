# SpaBench
SpaBench is a comprehensive spatial transcriptomics benchmarking workflow which implements 23 distinct methods (contain non-spatial and spatially-aware ones) and includes 360 different combinations of workflows. SpaBench begins from raw spatially-resolved gene count matrix, followed by normalization, feature selection, dimensionality reduction, and clustering. As a end point, SpaBench utilizes SingleR for cell type / spatial domain annotation and compute QC metrics, including both non-spatial ones such as NMI, HOM, COM, CHAOS, PAS, ASW, and spatially-aware IISC (Integrated Image–Spatial Cohesion).

## To run the whole workflow
```
module load snakemake
snakemake -j 36
```

## To run specific branch of workflow
```
snakemake -j60 -p -n   /path/todata/{CLR,LogNormalize,RC,SCTransform}/{vst,mvp,disp,markvariogram,moransi,nnSVG,SPARK,SPARK-X}/SpaGCN/06_SpaGCN_clustered.rds
```

