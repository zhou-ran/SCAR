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

  reference <- Matrix::as.matrix(reference[gene_intersect, ])


  cluster_inf <- dat@meta.data[,group.by]

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
