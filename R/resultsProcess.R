#' Processing peeling results
#'
#'
#' @param peel.results peeling results
#' @param posDT a data frame containing gene annotation information; a list
#'   component created by \code{dataPrep}.
#'
#' @return processed peeling results with a list of genes corresponding to each
#'   peeled region
#'
#' @details The \code{\link{peelingOneIterate}} function identifies (i) multiple loci across the genome
#'   where copy number gains are losses, and (ii) regions around those loci that also exhibit copy number
#'   changes.  Similarly, \code{\link{peelingTwoIterate}} identifies loci and surrounding regions that
#'   harbor copy number differences.  The output of either \code{\link{peelingOneIterate}} or
#'   \code{\link{peelingTwoIterate}} is processed to produce a tab-delimited text file that provides
#'   a summary of the peeling results.  Each column corresponds to a peeled locus, and the column contains
#'   the genomic location of the locus, the start and end positions of the surrounding peeled region,
#'   the mean copy number value (one matrix X) or difference of mean copy number values (two matrices X and Y),
#'   the cyclic shift-based p-value, and the names of the genes in the peeled region (in alphabetical order).
#'
#' @seealso \code{\link{dataPrep}}
#'
#' @import plyr
#' @export
resultsProcess = function(peel.results,posDT) {
  posDT$start = as.numeric(posDT$start)
  posDT$end = as.numeric(posDT$end)
  inner.fun = function(x) {
    x = as.numeric(x)
    chr.rows = which(posDT$chr == x[1])
    gene.rows1 = chr.rows[((posDT$start[chr.rows]-x[2])*(posDT$end[chr.rows]-x[3]) <= 0)]
    gene.rows2 = chr.rows[((posDT$start[chr.rows]-x[2])*(posDT$start[chr.rows]-x[3]) <= 0)]
    gene.rows3 = chr.rows[((posDT$end[chr.rows]-x[2])*(posDT$end[chr.rows]-x[3]) <= 0)]
    gene.rows = unique(c(gene.rows1,gene.rows2,gene.rows2))
    gene.result = data.frame(t(sort(rownames(posDT)[gene.rows])))
    colnames(gene.result) = paste0("G",1:length(gene.rows))
    return(gene.result)
  }
  gene.list = plyr::alply(peel.results[,c("chr","peelStart","peelStop")],1,inner.fun)
  result = t(cbind(peel.results,do.call(plyr::rbind.fill,gene.list)))
  return(result)
}
