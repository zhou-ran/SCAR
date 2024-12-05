#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj An SeuratData object.
#' @param idents the cluster id used for downstream analysis, Default: NULL
#' @param spots the spot id used for downstream analysis, Default: NULL
#' @param assay_1 the assay_1 used for downstream analysis, Default: 'Spatial'
#' @param gene_1 the gene_1 used for downstream analysis, Default: NULL
#' @param assay_2 the assay_2 used for downstream analysis, Default: 'Spatial'
#' @param gene_2 the gene_2 used for downstream analysis, Default: NULL
#' @param slot the slot used for downstream analysis, Default: 'data'
#' @param pt.size.factor pt.size.factor, Default: 1
#' @param max.cutoff The max cutoff for visualization, Default: NA
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
#'  \code{\link[checkmate]{checkFALSE}}
#'  \code{\link[Seurat]{reexports}}, \code{\link[Seurat]{SpatialPlot}}
#'  \code{\link[glue]{glue}}
#' @rdname SpatialInteractionPlot
#' @export
#' @importFrom checkmate assert_false
#' @importFrom Seurat Idents GetAssayData SpatialFeaturePlot
#' @importFrom glue glue
#'
#'
SpatialInteractionPlot <- function(obj,
                                   idents=NULL,
                                   spots=NULL,
                                   assay_1 = 'Spatial',
                                   gene_1 = NULL,
                                   assay_2 = 'Spatial',
                                   gene_2 = NULL,
                                   slot = 'data',
                                   pt.size.factor = 1,
                                   max.cutoff = NA,
                                   ...) {
  # only choose spots based on ident id or spot id
  checkmate::assert_false((!is.null(idents) & !is.null(spots)))
  if (!is.null(idents)) {
    obj <- obj[, Seurat::Idents(obj) %in% idents]
  }
  if (!is.null(spots)) {
    obj <- obj[, colnames(obj) %in% spots]
  }

  gene_label <- glue::glue('{gene_1}_{gene_2}')
  obj@meta.data[[gene_label]] <-
    Seurat::GetAssayData(obj, slot = slot, assay = assay_1)[gene_1, ] *
    Seurat::GetAssayData(obj, slot = slot, assay = assay_2)[gene_2, ]

  Seurat::SpatialFeaturePlot(
    obj,
    features = gene_label,
    crop = FALSE,
    stroke = NA,
    pt.size.factor = pt.size.factor,
    max.cutoff = max.cutoff
  ) -> p

  return(p)

}
