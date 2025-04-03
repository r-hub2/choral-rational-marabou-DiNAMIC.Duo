#' DNA copy number data for lung adenocarcinoma.
#'
#' A subset of the DNA copy number data from the TCGA lung
#' lung adenocarcinoma cohort.
#'
#' @format A numeric matrix of quantitative DNA copy number values with 3475 rows and 65 columns.
#'   Rows are indexed by genes, and columns are indexed by samples.  The entries are on the log
#'   ratio scale and were produced by the GISTIC pipeline.
#'
#' @source \url{https://gdac.broadinstitute.org/}
#'
#' @details The original data set downloaded from \url{https://gdac.broadinstitute.org/}
#'   contained gene-level segmented DNA copy number values for 24776 genes and 516 samples
#'   that were produced by the GISTIC pipeline.  Because of the size limitations on data
#'   in R packages, a random subset of genes and samples was selected.  By construction,
#'   the number of genes and samples in the lung adenocarcinoma data is different than the
#'   number of genes and samples in the lung squamous cell carcinoma data.
"luadSubset"