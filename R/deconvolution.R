#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param reference PARAM_DESCRIPTION
#' @param doublet_mode PARAM_DESCRIPTION, Default: 'full'
#' @param assay PARAM_DESCRIPTION, Default: 'Spatial'
#' @param ncores PARAM_DESCRIPTION, Default: 1
#' @param ... PARAM_DESCRIPTION
#' @return a Seurat Object
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
#' @seealso
#'  \code{\link[spacexr]{SpatialRNA}}, \code{\link[spacexr]{create.RCTD}}
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[Matrix]{colSums}}
#' @rdname run_rctd
#' @export
#' @importFrom spacexr SpatialRNA create.RCTD run.RCTD
#' @importFrom Seurat GetTissueCoordinates GetAssayData
#' @importFrom Matrix colSums
#'
#'

run_rctd <- function(obj,
                     reference,
                     doublet_mode = 'full',
                     assay = 'Spatial',
                     ncores = 1,
                     ...) {
  # check the class of reference object
  checkmate::assert_class(reference, 'Reference')

  rctd_obj <-
    spacexr::SpatialRNA(
      Seurat::GetTissueCoordinates(obj),
      Seurat::GetAssayData(obj, assay = assay, slot = 'counts'),
      Matrix::colSums(Seurat::GetAssayData(obj, assay = assay, slot = 'counts'))
    )

  rctd_obj <-
    spacexr::create.RCTD(spatialRNA = rctd_obj,
                         reference = reference,
                         max_cores = ncores)
  rctd_obj <-
    spacexr::run.RCTD(rctd_obj, doublet_mode = doublet_mode)

  decov_weight <-
    data.frame(spacexr::normalize_weights(rctd_obj@results$weights))

  obj@misc$rtcd <- decov_weight
  return(obj)
}


#' @title Perform de-convolution by DWLS
#' @description Perform de-convolution by DWLS, more detailed information refer to \url{https://github.com/dtsoucas/DWLS}
#' @param obj PARAM_DESCRIPTION
#' @param reference PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: 'Spatial'
#' @param group.by PARAM_DESCRIPTION, Default: 'cluster'
#' @param cutoff PARAM_DESCRIPTION, Default: 2
#' @param n_cell PARAM_DESCRIPTION, Default: 50
#' @param ... PARAM_DESCRIPTION
#' @return a Seurat Object
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[Matrix]{character(0)}}
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[Giotto]{enrich_deconvolution}}, \code{\link[Giotto]{spot_deconvolution}}
#' @rdname run_DWLS
#' @export
#' @importFrom Matrix as.matrix
#' @importFrom Seurat GetAssayData
#' @importFrom Giotto activeFeatType
#'

run_DWLS <- function(obj,
                     reference,
                     assay = 'Spatial',
                     group.by = 'cluster',
                     cutoff = 2,
                     n_cell = 50,
                     ...) {
  checkmate::assert_function(Giotto::activeFeatType)
  spots_raw_mtx <-
    expm1(Matrix::as.matrix(Seurat::GetAssayData(obj, assay = assay, slot = 'data')))

  spots_norm_mtx <-
    Matrix::as.matrix(Seurat::GetAssayData(obj, assay = assay, slot = 'data'))

  gene_intersect <-
    intersect(rownames(spots_raw_mtx), rownames(reference))

  reference <- Matrix::as.matrix(reference[gene_intersect,])


  cluster_inf <- dat@meta.data[, group.by]

  enrich_spot_proportion = Giotto:::enrich_deconvolution(
    expr = spots_raw_mtx,
    log_expr = spots_norm_mtx,
    cluster_info = cluster_inf,
    ct_exp = reference,
    cutoff = cutoff
  )

  resolution = (1 / n_cell)
  binarize_proportion = ifelse(enrich_spot_proportion >= resolution, 1, 0)

  spot_proportion <- Giotto:::spot_deconvolution(
    expr = spots_raw_mtx,
    cluster_info = cluster_inf,
    ct_exp = reference,
    binary_matrix = binarize_proportion
  )
  obj@misc$dwls <- t(data.frame(spot_proportion))
  return(obj)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param reference PARAM_DESCRIPTION
#' @param sc_meta PARAM_DESCRIPTION
#' @param assay PARAM_DESCRIPTION, Default: 'Spatial'
#' @param minCountGene PARAM_DESCRIPTION, Default: 0
#' @param minCountSpot PARAM_DESCRIPTION, Default: 0
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[checkmate]{checkClass}}
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[assertable]{assert_colnames}}
#'  \code{\link[CARD]{createCARDObject}}, \code{\link[CARD]{CARD_deconvolution}}
#' @rdname run_CARD
#' @export
#' @importFrom checkmate assert_class
#' @importFrom Seurat GetTissueCoordinates GetAssayData
#' @importFrom assertable assert_colnames
#' @importFrom CARD createCARDObject CARD_deconvolution
#'
run_CARD <- function(obj,
                     reference,
                     sc_meta,
                     assay = 'Spatial',
                     minCountGene = 0,
                     minCountSpot = 0,
                     ...) {
  checkmate::assert_class(reference, 'dgCMatrix')

  spatial_location <- Seurat::GetTissueCoordinates(obj)
  colnames(spatial_location) <- c('x', 'y')

  assertable::assert_colnames(
    data = sc_meta,
    colnames = c('cellID', 'cellType', 'sampleInfo'),
    quiet = T
  )

  CARD_obj = CARD::createCARDObject(
    sc_count = reference,
    sc_meta = sc_meta,
    spatial_count = Seurat::GetAssayData(obj, slot = 'counts', assay = assay),
    spatial_location = spatial_location,
    ct.varname = "cellType",
    ct.select = unique(sc_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = minCountGene,
    minCountSpot = minCountSpot
  )

  CARD_obj = CARD::CARD_deconvolution(CARD_object = CARD_obj)

  obj@misc$CARD <- CARD_obj@Proportion_CARD

  return(obj)
}

#'
#' #' @title FUNCTION_TITLE
#' #' @description FUNCTION_DESCRIPTION
#' #' @param obj PARAM_DESCRIPTION
#' #' @param reference PARAM_DESCRIPTION
#' #' @param sc_meta PARAM_DESCRIPTION
#' #' @param assay PARAM_DESCRIPTION, Default: 'Spatial'
#' #' @param ... PARAM_DESCRIPTION
#' #' @return OUTPUT_DESCRIPTION
#' #' @details DETAILS
#' #' @examples
#' #' \dontrun{
#' #' if(interactive()){
#' #'  #EXAMPLE1
#' #'  }
#' #' }
#' #' @seealso
#' #'  \code{\link[checkmate]{checkClass}}, \code{\link[checkmate]{checkLogical}}
#' #'  \code{\link[assertable]{assert_colnames}}
#' #'  \code{\link[SingleCellExperiment]{SingleCellExperiment-class}}, \code{\link[SingleCellExperiment]{colLabels}}
#' #'  \code{\link[scuttle]{logNormCounts}}
#' #'  \code{\link[scran]{modelGeneVar}}, \code{\link[scran]{getTopHVGs}}, \code{\link[scran]{scoreMarkers}}
#' #'  \code{\link[S4Vectors]{DataFrame-class}}, \code{\link[S4Vectors]{S4VectorsOverview}}
#' #'  \code{\link[Seurat]{reexports}}
#' #'  \code{\link[SpatialExperiment]{SpatialExperiment-class}}, \code{\link[SpatialExperiment]{SpatialExperiment}}
#' #'  \code{\link[SPOTlight]{SPOTlight}}
#' #' @rdname run_spotlight
#' #' @export
#' #' @importFrom checkmate assert_class assert_logical
#' #' @importFrom assertable assert_colnames
#' #' @importFrom SingleCellExperiment SingleCellExperiment colLabels
#' #' @importFrom scuttle logNormCounts
#' #' @importFrom scran modelGeneVar getTopHVGs scoreMarkers
#' #' @importFrom S4Vectors DataFrame
#' #' @importFrom Seurat GetTissueCoordinates GetAssayData
#' #' @importFrom SpatialExperiment SpatialExperiment
#' #' @importFrom SPOTlight SPOTlight
#' #'
#' run_spotlight <-
#'   function(obj, reference, sc_meta, assay = 'Spatial', ...) {
#'     checkmate::assert_class(reference, 'dgCMatrix')
#'     # checkmate::assert_logical(identical())
#'
#'     assertable::assert_colnames(
#'       data = sc_meta,
#'       colnames = c('cellID', 'cellType', 'sampleInfo'),
#'       quiet = T
#'     )
#'
#'     checkmate::assert_logical(identical(colnames(reference), sc_meta$cellID))
#'
#'     sce <-
#'       SingleCellExperiment::SingleCellExperiment(list(counts = reference))
#'
#'     sce$cell_type <- sc_meta$cellType
#'
#'     SingleCellExperiment::colLabels(sce) <- sce$cell_type
#'
#'     sce <- scuttle::logNormCounts(sce)
#'
#'     dec <- scran::modelGeneVar(sce)
#'     hvg <- scran::getTopHVGs(dec, n = 3000)
#'
#'     genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
#'     suppressMessages(require(S4Vectors))
#'     mgs <-
#'       scran::scoreMarkers(x=sce, subset.row = genes)
#'
#'     mgs_fil <- lapply(names(mgs), function(i) {
#'       x <- mgs[[i]]
#'       # Filter and keep relevant marker genes, those with AUC > 0.8
#'       x <- x[x$mean.AUC > 0.8, ]
#'       # Sort the genes from highest to lowest weight
#'       x <- x[order(x$mean.AUC, decreasing = TRUE), ]
#'       # Add gene and cluster id to the dataframe
#'       x$gene <- rownames(x)
#'       x$cluster <- i
#'       data.frame(x)
#'     })
#'
#'     mgs_df <- do.call(rbind, mgs_fil)
#'
#'     cd <- S4Vectors::DataFrame(Seurat::GetTissueCoordinates(obj))
#'
#'     spe <- SpatialExperiment::SpatialExperiment(
#'       assay = list(counts = as.matrix(
#'         Seurat::GetAssayData(obj, slot = 'count', assay = assay)
#'       )),
#'       colData = cd,
#'       spatialCoordsNames = c("imagerow", "imagecol")
#'     )
#'
#'     res <- SPOTlight::SPOTlight(
#'       x = sce,
#'       y = spe,
#'       groups = sce$cell_type,
#'       mgs = mgs_df,
#'       hvg = hvg,
#'       weight_id = "mean.AUC",
#'       group_id = "cluster",
#'       gene_id = "gene"
#'     )
#'
#'     obj@misc$SPOTlight <- as.data.frame(res$mat)
#'
#'     return(obj)
#'
#'   }
