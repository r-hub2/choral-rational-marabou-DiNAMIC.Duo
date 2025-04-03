#' Prepped data for DiNAMIC.Duo.
#'
#' A list produced by the dataPrep() function.
#'
#' @format A list containing three components for a common set of genes.
#'
#' @details The components of the list are named X, Y, and posDT.  They contain the DNA 
#'   copy number data for lung adenocarcinoma, the DNA copy number data for lung squamous 
#'   cell carcinoma, and gene position data, respectively.  This data object was produced 
#'   by applying \code{\link{dataPrep}} using the \code{\link{luadSubset}} and \code{\link{luscSubset}} 
#'   data sets.  This reduces run time when the package is compiled by CRAN, thus eliminating run 
#'   time errors.
"pD"