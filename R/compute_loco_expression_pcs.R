#!/usr/bin/env Rscript

library(argparse)
library(dplyr)

main <- function(args) {
  covar <- data.table::fread(args$covar, data.table = F, select = c(args$covar_sample_id, args$covariates))
  if (!is.null(args$covar_sample_id)) {
    covar <- dplyr::rename(covar, IID = !!as.symbol(args$covar_sample_id))
  }

  # gs://ea5_nea3/FinnGen_pseudobulk_all_cells.bed.gz
  expr_bed <- data.table::fread(args$expression_bed, data.table = F)
  expr_mat <- dplyr::select(expr_bed, -(chr:gene_id))

  # take samples that have covariates (non-PCA outlier)
  psam <- data.table::fread(args$psam, data.table = F)[, 1]
  samples <- intersect(colnames(expr_mat), psam)
  samples <- intersect(samples, covar$IID)
  samples <- samples[order(match(samples, colnames(expr_mat)))]
  covar <- dplyr::filter(covar, IID %in% samples)
  covar_mat <- model.matrix(~ . - 1, data = dplyr::select(covar, tidyselect::all_of(args$covariates)))
  n_samples <- length(samples)

  if (args$K > 0) {
    expr_mat <- as.matrix(expr_mat[, samples, drop = FALSE]) %>% t()

    if (args$decomposition == "FWL") {
      # regress out global covariates from expr_mat using FWL
      hat <- covar_mat %*% solve(t(covar_mat) %*% covar_mat) %*% t(covar_mat)
      I <- diag(nrow(covar_mat))
      expr_mat_res <- (I - hat) %*% expr_mat
    } else if (args$decomposition == "QR") {
      # regress out global covariates from expr_mat using QR
      qr_X <- qr(covar_mat)
      expr_mat_res <- qr.resid(qr_X, expr_mat)
    } else {
      stop("Unknown decomposition method")
    }
    rownames(expr_mat_res) <- samples

    if (!is.null(args$chromosome)) {
      loco_idx <- expr_bed$chr != args$chromosome
    } else {
      loco_idx <- rep(TRUE, nrow(expr_bed))
    }
    expr_pcs <- prcomp(expr_mat_res[, loco_idx], center = TRUE, scale. = TRUE, rank. = args$K)$x
    colnames(expr_pcs) <- paste0("GEX_PC", seq_len(ncol(expr_pcs)))
  } else {
    expr_pcs <- matrix(0, nrow = n_samples, ncol = 0, dimnames = list(samples, character(0)))
  }

  out_loco <-
    tibble::tibble(IID = samples) %>%
    dplyr::left_join(dplyr::select(covar, IID, dplyr::all_of(args$covariates)), by = "IID") %>%
    dplyr::bind_cols(expr_pcs)

  if (!is.null(args$chromosome)) {
    out_fname <- sprintf("%s.loco_%s.gex_pcs_covariates.tsv.gz", args$prefix, args$chromosome)
  } else {
    out_fname <- sprintf("%s.gex_pcs_covariates.tsv.gz", args$prefix)
  }
  data.table::fwrite(out_loco, out_fname, row.names = F, quote = F, na = "NA", sep = "\t")
}

parser <- ArgumentParser()

parser$add_argument("--prefix", type = "character", required = TRUE)
parser$add_argument("--expression-bed", type = "character", required = TRUE)
parser$add_argument("--psam", type = "character", required = TRUE)
parser$add_argument("--covar", type = "character", required = TRUE)
parser$add_argument("--covar-sample-id", type = "character", default = "IID")
parser$add_argument("--covariates", type = "character", nargs = "+", required = TRUE)
parser$add_argument("--chromosome", type = "character")
parser$add_argument("--decomposition", type = "character", choices = c("FWL", "QR"), default = "QR")
parser$add_argument("--K", "-K", type = "integer", default = 20)
parser$add_argument("--log", type = "character")

args <- parser$parse_args()

if (is.null(args$log)) {
  args$log <- paste0(args$prefix, ".log")
}

logfile <- file(args$log, open = "w")
sink(logfile, type = "output", split = TRUE)
sink(logfile, type = "message")
print(args)

print("Analysis started")
tryCatch(
  {
    main(args)
    print("Finished!")
  },
  error = function(e) {
    sink()
    message(as.character(e))
    sink(type = "message")
    stop(e)
  }
)
