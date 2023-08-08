#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj PARAM_DESCRIPTION
#' @param lr_network PARAM_DESCRIPTION
#' @param threshold_gene_exp PARAM_DESCRIPTION, Default: 0
#' @param slot_1 PARAM_DESCRIPTION, Default: 'counts'
#' @param assay_1 PARAM_DESCRIPTION, Default: 'Spatial'
#' @param slot_2 PARAM_DESCRIPTION, Default: 'counts'
#' @param assay_2 PARAM_DESCRIPTION, Default: 'Spatial'
#' @param prob PARAM_DESCRIPTION, Default: 'comb'
#' @param cores PARAM_DESCRIPTION, Default: 4
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
#'  \code{\link[stringr]{case}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[parallel]{mclapply}}
#'  \code{\link[cooccur]{cooccur}}
#'  \code{\link[plyr]{mapvalues}}
#' @rdname run_cooccur
#' @export
#' @importFrom stringr str_to_title
#' @importFrom tibble tibble
#' @importFrom parallel mclapply
#' @importFrom cooccur cooccur
#' @importFrom plyr mapvalues
#'
#'
run_cooccur <-
  function(obj,
           lr_network,
           threshold_gene_exp = 0,
           slot_1 = 'counts',
           assay_1 = 'Spatial',
           slot_2 = 'counts',
           assay_2 = 'Spatial',
           prob = 'comb',
           cores = 4,
           ...) {
    # first to check the dataframe format is from nichnet?
    # download from https://zenodo.org/record/3260758/files/lr_network.rds?download=1

    # replace human gene id into mouse gene id.

    lr_netword$from_mm <- stringr::str_to_title(lr_netword$from)
    lr_netword$to_mm <- stringr::str_to_title(lr_netword$to)

    # prepare matrix information
    hg_mtx <- Seurat::GetAssayData(obj,
                           slot = slot_1,
                           assay = assay_1)
    hg_binary_mtx <- +(hg_mtx > threshold_gene_exp)

    hg_gene_use <-
      rownames(hg_binary_mtx)[Matrix::rowSums(hg_binary_mtx) > 0]
    hg_gene_use <-
      hg_gene_use[hg_gene_use %in% unique(c(lr_netword$from, lr_netword$to))]


    mm_mtx <- Seurat::GetAssayData(dat, slot = slot_2, assay = assay_2)
    mm_binary_mtx <- +(mm_mtx > threshold_gene_exp)

    mm_gene_use <-
      rownames(mm_binary_mtx)[Matrix::rowSums(mm_binary_mtx) > 0]
    mm_gene_use <-
      mm_gene_use[mm_gene_use %in% unique(c(lr_netword$from_mm, lr_netword$to_mm))]


    final_pair_mtx <- tibble::tibble(
      from = c(lr_netword$from, lr_netword$from_mm),
      to = c(lr_netword$to_mm, lr_netword$to),
      source = c(lr_netword$source, lr_netword$source),
      database = c(lr_netword$database, lr_netword$database)
    )

    gene_all <- c(hg_gene_use, mm_gene_use)

    final_pair_mtx <-
      final_pair_mtx[final_pair_mtx$from %in% gene_all &
                       final_pair_mtx$to %in% gene_all, ]

    target_gene <- unique(c(final_pair_mtx$from, final_pair_mtx$to))
    final_mtx <-
      as.matrix(rbind(hg_binary_mtx[rownames(hg_binary_mtx) %in% target_gene,],
                      mm_binary_mtx[rownames(mm_binary_mtx) %in% target_gene,]))

    # split the gene into list
    index_spe_run <-
      split(final_pair_mtx, (seq_len(nrow(final_pair_mtx)) - 1) %/% 50)

    parallel::mclapply(index_spe_run, function(x) {
      x$label <- paste(x$from, x$to, sep = '_')
      tmp_mtx <- final_mtx[c(x$from, x$to), ]
      res_tmp <-
        cooccur::cooccur(
          mat = +tmp_mtx,
          prob = prob,
          spp_names = TRUE,
          thresh = FALSE
        )
      res_tmp <- res_tmp$results[, c(10, 11, 3:9)]
      res_tmp$label <-
        paste(res_tmp$sp1_name, res_tmp$sp2_name, sep = '_')
      res_tmp <- unique(res_tmp[res_tmp$label %in% x$label, ])
      res_tmp
    }, mc.cores = cores) -> res_para_comb



    final_pair_mtx$label <-
      paste(final_pair_mtx$from, final_pair_mtx$to, sep = '_')
    res_para_comb <- unique(do.call(rbind, res_para_comb))

    res_para_comb$sig <-
      ifelse(res_para_comb$p_lt >= 0.05, 0, -1) + ifelse(res_para_comb$p_gt >= 0.05, 0, 1)
    res_para_comb$pct_sp1 <-
      res_para_comb$obs_cooccur / res_para_comb$sp1_inc
    res_para_comb$pct_sp2 <-
      res_para_comb$obs_cooccur / res_para_comb$sp2_inc

    res_para_comb$source <-
      plyr::mapvalues(
        from = final_pair_mtx$label,
        to = final_pair_mtx$source,
        x = res_para_comb$label,
        warn_missing = FALSE
      )

    res_para_comb$database <-
      plyr::mapvalues(
        from = final_pair_mtx$label,
        to = final_pair_mtx$database,
        x = res_para_comb$label,
        warn_missing = FALSE
      )

    res_para_comb$coocer_rate <-
      apply(res_para_comb, 1, function(x) {
        min(as.numeric(x[['obs_cooccur']]) / as.numeric(x[['sp1_inc']]),
            as.numeric(x[['obs_cooccur']]) / as.numeric(x[['sp2_inc']]))
      })

    return(res_para_comb)
  }
