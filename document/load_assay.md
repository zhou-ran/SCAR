## load data
```r
library(SCAR)

dat <-
  load_data(
    data.dir = 'spatial_exp/hg',
    slice = 'HF1',
    assay_1 = 'hg',
    assay_2 = 'mm',
    assay_2_dir = 'spatial_exp/mm/'
  )
```

## Preprocessing datasets
```r
dat <-
  SCTransform(
    dat,
    assay = "hg",
    new.assay.name = "SCThg",
    return.only.var.genes = FALSE,
    verbose = FALSE
  )
dat <-
  SCTransform(
    dat,
    assay = "mm",
    new.assay.name = "SCTmm",
    return.only.var.genes = FALSE,
    verbose = FALSE
  )

dat <-
  RunPCA(
    dat,
    assay = "SCTmm",
    verbose = FALSE,
    reduction.name = 'mmpca',
    reduction.key = 'mmPC_'
  )

dat <- FindNeighbors(dat, reduction = "mmpca", dims = 1:30)
dat <- FindClusters(dat, verbose = FALSE, graph.name = 'SCTmm_snn')

dat <-
  RunUMAP(
    dat,
    reduction = "mmpca",
    dims = 1:30,
    reduction.name = 'mmumap',
    reduction.key = 'mmUMAP_'
  )

dat <-
  RunPCA(
    dat,
    assay = "SCThg",
    verbose = FALSE,
    reduction.name = 'hgpca',
    reduction.key = 'hgPC_'
  )

dat <- FindNeighbors(dat, reduction = "hgpca", dims = 1:30)
dat <- FindClusters(dat, verbose = FALSE, graph.name = 'SCThg_snn')
dat <-
  RunUMAP(
    dat,
    reduction = "hgpca",
    dims = 1:30,
    reduction.name = 'hgumap',
    reduction.key = 'hgUMAP_'
  )

saveRDS(dat, file = 'scar.Rds')

```
