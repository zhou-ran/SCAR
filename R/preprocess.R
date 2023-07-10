#' @title Title ST data preprocess
#' @description read 10X ST raw data and quality control
#' @param data.dir Directory containing the H5 file specified by filename and the image data in a subdirectory called spatial
#' @param filename Name of H5 file containing the feature barcode matrix
#' @param slice Name for the stored image of the tissue slice
#' @param assay_1 Name of the first assay
#' @param assay_2 Name of the second assay
#' @param assay_2_dir Directory containing the 10X market format expression matrix
#' @param ... PARAM_DESCRIPTION
#'
#' @return A Seurat object
#'
#' @details DETAILS
#' @examples
#' \dontrun{
#' dat <-
#' load_data(
#' data.dir = '/mnt/raid61/TNP_spatial/analysis/alignment/2021-1014-split/data/hg38/HF1/HF1/outs',
#' slice = 'HF1',
#' assay_1 = 'hg',
#' assay_2 = 'mm',
#' assay_2_dir = '/mnt/raid61/TNP_spatial/analysis/alignment/2021-1014-split/data/mm10/HF1/HF1/outs/filtered_feature_bc_matrix/'
#' )
#' }
#' @seealso
#'  \code{\link[Seurat]{Load10X_Spatial}}, \code{\link[Seurat]{Read10X}}, \code{\link[Seurat]{reexports}}
#'  \code{\link[checkmate]{checkTRUE}}
#'  \code{\link[glue]{glue}}
#' @rdname load_data
#' @export
#' @importFrom Seurat Load10X_Spatial Read10X CreateAssayObject
#' @importFrom checkmate assert_true
#' @importFrom glue glue
#'
#'
load_data <-
  function(data.dir,
           filename = 'filtered_feature_bc_matrix.h5',
           slice = 'slice1',
           assay_1 = 'hg',
           assay_2 = NULL,
           assay_2_dir = NULL,
           ...) {
    message(paste(Sys.time(), "[INFO] load the dataset"))
    dat <-
      Seurat::Load10X_Spatial(data.dir = data.dir,
                              slice = slice,
                              assay = assay_1)

    checkmate::assert_true(is.null(assay_2) == is.null(assay_2_dir))


    if (!is.null(assay_2) & !is.null(assay_2_dir)) {
      message(paste(
        Sys.time(),
        glue::glue("[INFO] add {assay_2_dir} into assay ('{assay_2}')")
      ))
      assay_2_obj <-
        Seurat::Read10X(data.dir =  assay_2_dir)

      dat[[assay_2]] <- Seurat::CreateAssayObject(assay_2_obj)

    }


    return(dat)

  }
