# Deconvolution


## deconvolution using RCTD
```r
## reference: https://github.com/dmcable/spacexr
ref <-
  readRDS('spacexr_2000ds.Rds')

## run RCTD
dat <- TBD::run_rctd(
  obj = dat,
  reference = ref,
  doublet_mode = 'full',
  ncores = 24,
  assay = 'mm'
)

```


## DWLS deconvolution
```r
### prepare signature matrix
seu <- CreateSeuratObject(counts = ref@counts)
seu$cell_type <- ref@cell_types[colnames(seu)]
seu <- Seurat::NormalizeData(seu)
Idents(seu) <- 'cell_type'
de_gene <- FindAllMarkers(seu, only.pos = TRUE, test.use = "bimod")

pval.cutoff <- 0.01
diff.cutoff <- 0.5

DEGenes <-
  de_gene$gene[intersect(which(de_gene$p_val_adj < pval.cutoff),
                         which(de_gene$avg_log2FC > diff.cutoff))]
nonMir = grep("MIR|Mir", DEGenes, invert = T)
de_gene_use <- de_gene[which(de_gene$gene %in% DEGenes[nonMir]),]
norm_mtx <- expm1(Seurat::GetAssayData(seu, 'data'))

lapply(levels(Idents(seu)), function(x) {
  Matrix::rowMeans(norm_mtx[DEGenes, Idents(seu) == x])
}) -> sig_mtx
sig_mtx <- do.call(cbind, sig_mtx)
colnames(sig_mtx) <- levels(Idents(seu))

dat <-
  run_DWLS(
    obj = dat,
    reference = sig_mtx,
    assay = 'SCTmm',
    group.by = 'SCTmm_snn_res.0.8'
  )
```

## CARD testing
```r
sc_meta <-
  data.frame(
    cellID = colnames(ref@counts),
    cellType = ref@cell_types,
    sampleInfo = 'sample1'
  )

dat <- run_CARD(
  obj = dat,
  reference = Matrix::Matrix(ref@counts, sparse = TRUE),
  sc_meta = sc_meta,
  assay = 'SCTmm'
)
```

## SPOTlight
```r
dat <-
  run_spotlight(
    obj = dat,
    reference = ref@counts,
    sc_meta = sc_meta,
    assay = 'SCTmm'
  )

### cooccur interaction analysis
res_para_comb <- run_cooccur(
  dat,
  lr_network = lr_netword,
  assay_1 = 'SCThg',
  assay_2 = 'SCTmm',
  cores = 24
)
```