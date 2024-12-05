#' @title calculate geneSet score
#' @description Calculate gene set score based on multiple methods
#' @param obj The SeuratData object for analysis.
#' @param assay The assasy used for score analysis.
#' @param signatures GeneSet used which could be a GeneSetCollection object, name of Gmt used in GSEABase and one of c('KEGG_metabolism_nc', 'REACTOME_metabolism'), Default: 'KEGG'
#' @param costume_label The label of score, Default: 'costume'
#' @param ncores The core used for downstream analysis, Default: 1
#' @param method_use Method used for calculating score, including c("AUC","VISION","ssgsea"), Default: 'AUC'
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
#'  \code{\link[checkmate]{checkChoice}}, \code{\link[checkmate]{checkClass}}
#'  \code{\link[GSEABase]{import/export}}
#'  \code{\link[tools]{fileutils}}
#'  \code{\link[AUCell]{AUCell_buildRankings}}, \code{\link[AUCell]{AUCell_calcAUC}}, \code{\link[AUCell]{aucellResults-class}}
#'  \code{\link[glue]{glue}}
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[VISION]{c("Vision", "VISION", "VISION")}}, \code{\link[VISION]{analyze,Vision-method}}
#' @rdname calculate_score
#' @export
#' @importFrom checkmate check_choice check_class
#' @importFrom GSEABase getGmt
#' @importFrom tools file_path_sans_ext
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @importFrom glue glue
#' @importFrom Matrix colSums
#' @importFrom VISION Vision analyze
#'

calculate_score <- function(obj,
                            assay,
                            signatures = 'KEGG',
                            costume_label = 'costume',
                            ncores = 1,
                            method_use = 'AUC',
                            ...) {
  #check the name of assay
  assay_labels <- names(obj@assays)
  checkmate::check_choice(assay, assay_labels)

  # check the class of reference object
  if (class(signatures) == 'character') {
    if (file.exists(signatures, showWarnings = FALSE)) {
      geneSets <- GSEABase::getGmt(signatures)
      signatures <- tools::file_path_sans_ext(basename(signatures))
      gmtfile <- signatures
    } else {
      checkmate::check_choice(
        x = signatures,
        choices = c('KEGG_metabolism_nc', 'REACTOME_metabolism')
      )

      gmtfile <-
        system.file("data", "{signatures}.gmt", package = packageName(environment(load_data)))
      geneSets <- GSEABase::getGmt(gmtfile)
    }

  } else {
    checkmate::check_class(signatures, 'GeneSetCollection')
    gmtfile <- NULL
    geneSets <- signatures
    signatures <- costume_label
  }

  # load gene signatures based on AUC
  if (method_use == 'AUC') {
    cells_rankings <-
      AUCell::AUCell_buildRankings(
        GetAssayData(obj, slot = 'count', assay = assay),
        splitByBlocks = TRUE,
        plotStats = FALSE
      )

    cells_AUC <-
      AUCell::AUCell_calcAUC(geneSets,
                             cells_rankings,
                             nCores = ncores,
                             verbose = FALSE)

    signature_exp <- data.frame(AUCell::getAUC(cells_AUC))
    colnames(signature_exp) <-
      colnames(GetAssayData(obj, slot = 'count', assay = assay))

    obj[[glue::glue('score{signatures}')]] <-
      CreateAssayObject(data = as.matrix(signature_exp))
  }

  # load gene signatures based on VISION
  if (method_use == "VISION") {
    n.umi <-
      Matrix::colSums(GetAssayData(obj, slot = 'count', assay = assay))

    vis <-
      VISION::Vision(t(t(
        GetAssayData(obj, slot = 'count', assay = assay)
      ) / n.umi) * median(n.umi), signatures = gmtFile)

    vis <- VISION::analyze(vis)

    signature_exp <- data.frame(t(vis@SigScores))
    colnames(signature_exp) <-
      colnames(GetAssayData(obj, slot = 'count', assay = assay))

    obj[[glue::glue('score{costume_label}')]] <-
      CreateAssayObject(data = as.matrix(signature_exp))
  }

  # load gene signatures based on ssgsea
  if (method_use == "ssgsea") {
    gsva_es <-
      gsva(
        as.matrix(obj@assays$SCThg@counts),
        geneSets,
        method = c("ssgsea"),
        kcdf = c("Poisson"),
        parallel.sz = 2
      ) #
    signature_exp <- data.frame(gsva_es)
    colnames(signature_exp) <-
      colnames(GetAssayData(obj, slot = 'count', assay = assay))

    obj[[glue::glue('score{signatures}')]] <-
      CreateAssayObject(data = as.matrix(signature_exp))

  }
  return(obj)
}
