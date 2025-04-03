#' A Function to Apply the Peeling Algorithm for Two Copy Number Matrices
#'
#' This function iteratively applies the peelingTwo function, thereby identifying multiple
#'
#' differences across the genome between a two cohorts.  Gains and losses should be analyzed separately.
#'
#' @param X A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
#'
#' @param Y A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
#'
#' @param posDT A data frame containing genomic position information for the genes in X.
#'
#' @param gain A logical value indicating whether gains (TRUE) or losses (FALSE) will be peeled.
#'
#' Default = TRUE.
#'
#' @param nullDist An empirical null distribution produced by the cyclic shift algorithm.
#'
#' Default = NULL.
#'
#' @param threshold A tuning parameter that controls the size of the peeled region.  Rows of X
#'
#' and Y with mean copy number differences less than threshold will not be peeled.  Default = NULL.
#'
#' @param numIters The number of times peelingTwo will be iterated.  Default = 5.
#'
#' @return A list containing two elements: X, Y, and interval.  X and Y are updated versions of the input
#' copy number matrices in which the peak difference at k has been removed, and interval is genomic region
#' containing k.  By construction, interval cannot extend beyond the chromosome arm containing k.
#'
#' @details The \code{\link{peelingTwo}} function applies the peeling procedure for two cohorts to "nullify"
#'   entries in two copy number matrices X and Y that give rise to the most significant copy number
#'   difference.  Because tumor genomes may contain multiple regions that harbor copy number differences,
#'   it is important to apply the peeling procedure for two cohorts iteratively, thereby identifying
#'   multiple markers and surrounding genomic regions.
#' 
#' @examples luad=pD[["X"]]
#'
#'  lusc=pD[["Y"]]
#'
#'  posDT=pD[["posDT"]]
#'
#'  gain = TRUE
#'
#'  nullDist = NULL
#'
#'  threshold = NULL
#'
#'  numIters = 3
#'
#'  out=peelingTwoIterate(X=luad,Y=lusc,posDT=posDT,gain=TRUE,nullDist=NULL,threshold=NULL,numIters=3)   
#'
#' @export

peelingTwoIterate = function (X, Y, posDT, gain = TRUE, nullDist = NULL, threshold = NULL, numIters = 5)
	{
	output = data.frame(chr = rep(NA, numIters), peelStart = rep(NA, numIters), peelStop = rep(NA, numIters),
		peakLoc = rep(NA, numIters), peakVal = rep(NA, numIters), peakPVal = rep(NA, numIters))
	gainLossInd = (-1)^(1 + gain)
	X1 = X * gainLossInd
	Y1 = Y * gainLossInd
	for (i in 1:numIters)
		{
		tempK = which.max(rowMeans(X1) - rowMeans(Y1))
		tempPeakVal = gainLossInd * signif((rowMeans(X1) - rowMeans(Y1))[tempK], 4)
		tempPeakPVal = ifelse(is.null(nullDist), NA, (1 + sum(nullDist > tempPeakVal))/length(nullDist) * (gainLossInd == 1) +
			(1 + sum(nullDist < tempPeakVal))/length(nullDist) * (gainLossInd == -1))
		tempPeakPval = min(tempPeakPVal, 1)
		tempPeel = peelingTwo(X = X1, Y = Y1, posDT = posDT, k = tempK, threshold = threshold)
		output[i,] = c(posDT[tempK, 1], tempPeel[[3]][1], tempPeel[[3]][2], posDT[tempK, 2], 
			tempPeakVal, tempPeakPVal)
		X1 = tempPeel[[1]]
		Y1 = tempPeel[[2]]
		}
	return(output)
	}
