#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj A SeuratData object.
#' @param lr_network LR database which was provided by SCAR
#' @param idents Only check the given cluster which should be one of cluster in `Idents(obj)`, Default: NULL
#' @param spots Only check the given spots, Default: NULL
#' @param threshold_gene_exp The threshold for considering a gene as expressed, Default: 0
#' @param slot_1 the slot of dataset1 used for LR analysis, Default: 'data'
#' @param assay_1 the assay of dataset1 used for LR analysis, Default: 'Spatial'
#' @param slot_2 the slot of dataset2 used for LR analysis, Default: 'data'
#' @param assay_2 the assay of dataset2 used for LR analysis, Default: 'Spatial'
#' @param cores the core for analysis, Default: 4
#' @param prob Should co-occurrence probabilities be calculated using the hypergeometric distribution (prob="hyper") or the combinatorics approach from Veech 2013 (prob="comb"). Default: 'comb'
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
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[stringr]{case}}
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[parallel]{mclapply}}
#'  \code{\link[cooccur]{cooccur}}
#'  \code{\link[plyr]{mapvalues}}
#' @rdname run_cooccur
#' @export
#' @importFrom checkmate assert_false
#' @importFrom Seurat Idents GetAssayData
#' @importFrom stringr str_to_title
#' @importFrom Matrix rowSums
#' @importFrom tibble tibble
#' @importFrom parallel mclapply
#' @importFrom cooccur cooccur
#' @importFrom plyr mapvalues
#'
run_cooccur <-
  function(obj,
           lr_network,
           idents=NULL,
           spots=NULL,
           threshold_gene_exp = 0,
           assay_1 = 'Spatial',
           slot_1 = 'counts',
           assay_2 = 'Spatial',
           slot_2 = 'counts',
           prob = 'comb',
           cores = 4,
           ...) {
    # only choose spots based on ident id or spot id
    checkmate::assert_false((!is.null(idents) & !is.null(spots)))
    if (!is.null(idents)) {
      obj <- obj[, Seurat::Idents(obj) %in% idents]
    }
    if (!is.null(spots)) {
      obj <- obj[, colnames(obj) %in% spots]
    }
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


    mm_mtx <- Seurat::GetAssayData(obj, slot = slot_2, assay = assay_2)
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

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param obj A SeuratData object.
#' @param lr_network LR database which was provided by SCAR
#' @param idents Only check the given cluster which should be one of cluster in `Idents(obj)`, Default: NULL
#' @param spots Only check the given spots, Default: NULL
#' @param threshold_gene_exp The threshold for considering a gene as expressed, Default: 0
#' @param slot_1 the slot of dataset1 used for LR analysis, Default: 'data'
#' @param assay_1 the assay of dataset1 used for LR analysis, Default: 'Spatial'
#' @param slot_2 the slot of dataset2 used for LR analysis, Default: 'data'
#' @param assay_2 the assay of dataset2 used for LR analysis, Default: 'Spatial'
#' @param cores the core for analysis, Default: 4
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
#'  \code{\link[Seurat]{reexports}}
#'  \code{\link[stringr]{case}}
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[tibble]{tibble}}
#'  \code{\link[parallel]{mclapply}}
#' @rdname run_iascore
#' @export
#' @importFrom checkmate assert_false
#' @importFrom Seurat Idents GetAssayData
#' @importFrom stringr str_to_title
#' @importFrom Matrix rowSums rowMeans
#' @importFrom tibble tibble
#' @importFrom parallel mclapply
#'
run_iascore <-
  function(obj,
           lr_network,
           idents=NULL,
           spots=NULL,
           threshold_gene_exp = 0,
           slot_1 = 'data',
           assay_1 = 'Spatial',
           slot_2 = 'data',
           assay_2 = 'Spatial',
           cores = 4,
           ...) {
    # only choose spots based on ident id or spot id
    checkmate::assert_false((!is.null(idents) & !is.null(spots)))
    if (!is.null(idents)) {
      obj <- obj[, Seurat::Idents(obj) %in% idents]
    }
    if (!is.null(spots)) {
      obj <- obj[, colnames(obj) %in% spots]
    }

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


    mm_mtx <-
      Seurat::GetAssayData(obj, slot = slot_2, assay = assay_2)
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
    gene_all <- c(hg_gene_use, mm_gene_use)

    final_pair_mtx <-
      final_pair_mtx[final_pair_mtx$from %in% gene_all &
                       final_pair_mtx$to %in% gene_all, ]


    mtx_use <-
      rbind(mm_mtx[rownames(mm_mtx) %in% c(final_pair_mtx$from, final_pair_mtx$to),],
            hg_mtx[rownames(hg_mtx) %in% c(final_pair_mtx$from, final_pair_mtx$to),])


    mean_exp <- Matrix::rowMeans(mtx_use)

    exp_rate <- Matrix::rowMeans(+(mtx_use > 0))

    interaction_df <- data.frame(
      from = final_pair_mtx$from,
      to = final_pair_mtx$to,
      source = paste(final_pair_mtx$database, final_pair_mtx$source, sep =
                       ':'),
      score = mean_exp[final_pair_mtx$from] * mean_exp[final_pair_mtx$to],
      rate = exp_rate[final_pair_mtx$from] * exp_rate[final_pair_mtx$to]
    )

    ## generate a random interaction pair to simulate a distribution

    simu_distri <- parallel::mclapply(1:10000, function(i) {
      gene_pair <- sample(rownames(mtx_use), 2)
      return(prod(Matrix::rowMeans(mtx_use[gene_pair,])))
    }, mc.cores = cores)

    simu_distri <- unlist(simu_distri)


    interaction_df$pval <-
      unlist(lapply(interaction_df$score, function(x) {
        1 - sum(+(x > simu_distri)) / (length(simu_distri) + 1)
      }))
    
    interaction_df$padj <- p.adjust(interaction_df$pval)
    interaction_df$label <- paste(interaction_df$from, interaction_df$to, sep = '_')

    return(interaction_df)
  }
