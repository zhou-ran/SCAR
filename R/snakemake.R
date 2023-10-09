
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dir PARAM_DESCRIPTION
#' @param config_list PARAM_DESCRIPTION, Default: NULL
#' @param overwrite PARAM_DESCRIPTION, Default: FALSE
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
    # config_lst <-       list(
    #   'output' = 'data',
    #   'samples' = c('LL'),
    #   'FASTQ' = 'fq',
    #   'STAR' = list('mix_index' = 'index'),
    #                 'cores' = as.integer(8)),
    #   'soft' = list(
    #     'STAR' = 'STAR',
    #     'seqkit' = 'seqkit',
    #     'split_script' = 'bin/split_read_id.py',
    #     'seqkt' = 'seqkt',
    #     'python' = 'python'
    #   )
    # )
    yaml::write_yaml(x = config_list,
                     file = file.path(dir, 'config.yaml'))
  }
