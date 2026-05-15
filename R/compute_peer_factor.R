#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(methods)

# ref: https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/src/run_PEER.R
# ref: https://github.com/immunogenomics/sceQTL/blob/main/scripts/pseudobulk/02_normalize.R
main <- function(args) {
  covar <- data.table::fread(args$covar, data.table = F, select = c(args$covar_sample_id, args$covariates)) %>%
    dplyr::rename(IID = !!as.symbol(args$covar_sample_id))

  # gs://ea5_nea3/FinnGen_pseudobulk_all_cells.bed.gz
  expr_bed <- data.table::fread(args$expression_bed, data.table = F)
  expr_mat <- dplyr::select(expr_bed, -(chr:gene_id)) %>%
    as.matrix()

  # take samples that have covariates (non-PCA outlier)
  psam <- data.table::fread(args$psam, data.table = F)[, 1]
  samples <- intersect(colnames(expr_mat), psam)
  samples <- intersect(samples, covar$IID)

  if (!is.null(args$subset_samples)) {
    subset_samples <- read.table(args$subset_samples, header = F, as.is = T)[, 1]
    samples <- intersect(samples, subset_samples)
  }
  covar <- dplyr::filter(covar, IID %in% samples)
  n_samples <- length(samples)

  # PEER normalization: regress out latent variables (Stegle 2010, PLoS Comp Bio)
  # Based on GTEx 2017 Nature:
  # https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl
  # The number of PEER factors was selected as function of sample size (N):
  # * 15 factors for N < 150
  # * 30 factors for 150 ≤ N < 250
  # * 45 factors for 250 ≤ N < 350
  # * 60 factors for N ≥ 350
  if (!is.null(args$K)) {
    K <- args$K
  } else if (n_samples < 150) {
    K <- 15
  } else if (n_samples < 250) {
    K <- 30
  } else if (n_samples < 350) {
    K <- 45
  } else if (n_samples >= 350) {
    K <- 60
  }

  if (K > 0) {
    model <- peer::PEER()
    peer::PEER_setPhenoMean(model, as.matrix(t(expr_mat)))
    # peer::PEER_setAdd_mean(model, TRUE)
    peer::PEER_setNk(model, K)
    peer::PEER_setCovariates(model, as.matrix(covar))
    peer::PEER_setNmax_iterations(model, args$max_iterations)
    peer::PEER_update(model)

    X <- peer::PEER_getX(model) # samples x PEER factors
    rownames(X) <- colnames(expr_mat)
    colnames(X) <- paste0("PEER", seq(ncol(X)))
    A <- peer::PEER_getAlpha(model) # PEER factors x 1
    colnames(A) <- "Alpha"
    R <- t(peer::PEER_getResiduals(model)) # genes x samples
    colnames(R) <- colnames(expr_mat)

    if (all(is.na(A))) {
      stop("Model doesn't converge.")
    }
  } else {
    print("No PEER factors are computed.")
    X <- matrix(0, nrow = n_samples, ncol = 0, dimnames = list(samples, character(0)))
    A <- data.frame(Alpha = NA)
    R <- NULL
  }

  out_peer <- tibble::as_tibble(X, rownames = "IID") %>%
    dplyr::inner_join(dplyr::select(covar, IID, dplyr::all_of(args$covariates)), by = "IID")

  out_alpha <- as.data.frame(A) %>%
    dplyr::mutate(Relevance = 1.0 / Alpha)

  out_bed <- dplyr::bind_cols(
    dplyr::select(expr_bed, dplyr::matches("chr$"), start, end, gene_id),
    R
  )

  data.table::fwrite(out_peer, sprintf("%s.peer_factors_covariates.tsv.gz", args$prefix), row.names = F, quote = F, na = "NA", sep = "\t")
  data.table::fwrite(out_alpha, sprintf("%s.peer_alpha.tsv.gz", args$prefix), row.names = F, quote = F, na = "NA", sep = "\t")
  data.table::fwrite(out_bed, sprintf("%s.peer_residuals.bed.gz", args$prefix), row.names = F, quote = F, na = "NA", sep = "\t")
}

parser <- ArgumentParser()

parser$add_argument("--prefix", type = "character", required = TRUE)
parser$add_argument("--expression-bed", type = "character", required = TRUE)
parser$add_argument("--psam", type = "character", required = TRUE)
parser$add_argument("--covar", type = "character", required = TRUE)
parser$add_argument("--covar-sample-id", type = "character", default = "IID")
parser$add_argument("--covariates", type = "character", nargs = "+", required = TRUE)
parser$add_argument("--subset-samples", type = "character")
parser$add_argument("--K", "-K", type = "integer")
parser$add_argument("--max-iterations", type = "integer", default = 10000)

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
