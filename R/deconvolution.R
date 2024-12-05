#' @title An integrative function for performing deconvolution analysis by RCTD.
#' @description An integrative function for performing deconvolution analysis by RCTD.
#' @param obj The SeuratData object for analysis.
#' @param reference An RCTD-created reference. 
#' @param if TRUE, uses RCTD doublet mode weights. Otherwise, uses RCTD full mode weights. Default: 'TRUE'
#' @param assay The assasy used for deconvolution analysis, Default: 'Spatial'
#' @param ncores The cores used for deconvolution analysis, Default: 1
#' @param ... PARAM_DESCRIPTION
#' @return a Seurat Object
#' @details An integrative function for performing deconvolution analysis by RCTD.
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
                     doublet_mode = TRUE,
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
#' @param obj The SeuratData object for analysis.
#' @param reference A raw count matrix which rownames is the gene id.
#' @param assay The assasy used for deconvolution analysis, Default: 'Spatial'
#' @param group.by The cluster id used for analysis, Default: 'cluster'
#' @param cutoff cutoff used in Giotto deconvolution, Default: 2
#' @param n_cell the number of spots expressed, Default: 50
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
#' @param obj The SeuratData object for analysis.
#' @param reference A dgCMatrix which rownames is the gene id.
#' @param sc_meta A data.frome whcih contains c('cellID', 'cellType', 'sampleInfo') colnames.
#' @param assay The assasy used for deconvolution analysis, Default: 'Spatial'
#' @param minCountGene minCountGene, Default: 0
#' @param minCountSpot minCountSpot, Default: 0
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

