#!/usr/bin/env Rscript

library(ArchR)
library(argparse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

read_narrowPeaks <- function(narrowPeakFile, summits = FALSE, extendSummits = 250, mergeSummits = FALSE, maxPeaks = NULL, q_threshold = NULL) {
    out <- data.table::fread(narrowPeakFile) %>%
        dplyr::mutate(
            start = if (summits) {
                V2 + V10 + 1
            } else {
                V2 + 1
            },
            end = if (summits) {
                V2 + V10 + 1
            } else {
                V3
            }
        )
    if (!is.null(maxPeaks)) {
        sorted_counts <- sort(out$V11, decreasing = TRUE)
        count_threshold <- ifelse(length(sorted_counts) > maxPeaks, sorted_counts[maxPeaks], sorted_counts[length(sorted_counts)])
        out <- dplyr::filter(out, V11 >= count_threshold)
    }
    if (!is.null(q_threshold)) {
        out <- dplyr::filter(out, V9 > -log10(q_threshold))
    }

    out <- GRanges(
        out$V1,
        IRanges(out$start, out$end),
        name = out$V4,
        score = out$V5,
        quantile_score = ArchR:::.getQuantiles(out$V5),
        nlog10p = out$V8,
        nlog10q = out$V9,
        relative_summit_pos = out$V10,
        count = out$V11
    )
    if (summits) {
        out <- GenomicRanges::resize(out, extendSummits * 2 + 1, "center")
    }

    out <- sort(sortSeqlevels(out))
    if (mergeSummits) {
        out <-
            bedtoolsr::bt.merge(out, c = "6,12", o = "distinct,max") %>%
            dplyr::rename(seqnames = V1, start = V2, end = V3, name = V4, count = V5) %>%
            GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    }
    return(out)
}

iterative_overlap_peak_merge <- function(
    peaks,
    score_column,
    promoterRegion = c(2000, 100),
    verbose = TRUE,
    threads = ArchR::getArchRThreads()) {
    ArchR::addArchRGenome("hg38")
    genomeAnnotation <- ArchR::getGenomeAnnotation()
    geneAnnotation <- ArchR::getGeneAnnotation()
    geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
    genomeAnnotation <- ArchR:::.validGenomeAnnotation(genomeAnnotation)
    geneAnnotation <- ArchR:::.validGeneAnnoByGenomeAnno(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation)

    logFile <- createLogFile("addReproduciblePeakSet")
    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)

    #####################################################
    # BSgenome for Add Nucleotide Frequencies!
    #####################################################
    ArchR:::.requirePackage("Biostrings", source = "bioc")
    BSgenome <- eval(parse(text = genomeAnnotation$genome))
    BSgenome <- validBSgenome(BSgenome)

    peaks <- sort(sortSeqlevels(peaks))
    peaks <- subsetByOverlaps(peaks, genomeAnnotation$blacklist, invert = TRUE)
    peaks <- subsetByOverlaps(peaks, genomeAnnotation$chromSizes, type = "within")
    peaks <- ArchR:::.fastAnnoPeaks(peaks, BSgenome = BSgenome, geneAnnotation = geneAnnotation, promoterRegion = promoterRegion, logFile = logFile)
    peaks <- peaks[which(mcols(peaks)$N < 0.001)] # Remove N Containing Peaks
    mcols(peaks)$N <- NULL # Remove N Column

    peaks <- ArchR::nonOverlappingGR(peaks, by = score_column, decreasing = TRUE)
    peaks <- sort(sortSeqlevels(peaks))
    return(peaks)
}

# peaks <- Sys.glob("*.narrowPeak") %>%
#     purrr::map(read_narrowPeaks, summits = TRUE, mergeSummits = TRUE, maxPeaks = 150000) %>%
#     GenomicRanges::GRangesList() %>%
#     unlist()
# merged_peaks <- iterative_overlap_peak_merge(peaks, "count")
#
# as.data.frame(merged_peaks) %>%
#     dplyr::select(seqnames, start, end, name, count, strand) %>%
#     data.table::fwrite("integrated_atac_batch1_4b.iterative_overlap_merged.bed", quote = F, col.names = F, na = "NA", sep = "\t")

main <- function(args) {
    peaks <- purrr::map(args$narrowPeak, read_narrowPeaks, summits = TRUE, mergeSummits = TRUE, maxPeaks = args$maxPeaks) %>%
        GenomicRanges::GRangesList() %>%
        unlist()
    merged_peaks <- iterative_overlap_peak_merge(peaks, "count")

    df_peaks <- as.data.frame(merged_peaks)
    data.table::fwrite(df_peaks, paste0(args$prefix, ".iterative_merged.tsv"), quote = F, na = "NA", sep = "\t")

    dplyr::select(df_peaks, seqnames, start, end, name, count, strand) %>%
        data.table::fwrite(paste0(args$prefix, ".iterative_merged.bed"), quote = F, col.names = F, na = "NA", sep = "\t")
}

parser <- ArgumentParser()
parser$add_argument("--prefix", type = "character", required = TRUE)
parser$add_argument("--narrowPeak", type = "character", nargs = "+", required = TRUE)
parser$add_argument("--maxPeaks", type = "integer", default = 150000)

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
