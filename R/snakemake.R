
#' @title Generate snakemake configure folder for demultiplexing reads
#' @description Generate snakemake configure folder for demultiplexing reads
#' @param dir the name of folder to save the snakemake configure.
#' @param config_list the configure ionformation to generate yaml file, Default: NULL
#' @param overwrite Whether to overwrite the folder if existed, Default: FALSE
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' config_lst <-       list(
#'   'output' = 'data',
#'   'samples' = c('LL'),
#'   'FASTQ' = 'fq',
#'   'STAR' = list('mix_index' = 'index',
#'                 'cores' = as.integer(8)),
#'   'soft' = list(
#'     'STAR' = 'STAR',
#'     'seqkit' = 'seqkit',
#'     'split_script' = 'bin/split_read_id.py',
#'     'seqkt' = 'seqkt',
#'     'python' = 'python'
#'   )
#' )
#'
#' init_smk_script(dir = 'smk_test/',
#'                 config_list = config_lst,
#'                 overwrite = T)
#' }
#' @seealso
#'  \code{\link[checkmate]{checkPathForOutput}}, \code{\link[checkmate]{checkList}}
#'  \code{\link[yaml]{write_yaml}}
#' @rdname init_smk_script
#' @export
#' @importFrom checkmate checkPathForOutput check_list
#' @importFrom yaml write_yaml
#'

init_smk_script <-
  function(dir,
           config_list = NULL,
           overwrite = FALSE,
           ...) {
    if (isFALSE(overwrite)) {
      checkmate::checkPathForOutput(dir)
    }

    checkmate::check_list(config_list)

    dir.create(file.path(dir, 'bin'),
               recursive = TRUE,
               showWarnings = FALSE)

    smk_file <-
      system.file("data", "smk.py", package = packageName(environment(load_data)))
    split_script <-
      system.file("data", "split_read_id.py", package = packageName(environment(load_data)))

    file.copy(from = smk_file,
              to = dir,
              overwrite = overwrite)

    file.copy(from = split_script,
              to = file.path(dir, 'bin'),
              overwrite = overwrite)

    yaml::write_yaml(x = config_list,
                     file = file.path(dir, 'config.yaml'))
  }
