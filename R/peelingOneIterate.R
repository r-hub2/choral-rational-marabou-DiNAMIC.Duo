#' A Function to Apply the Peeling Algorithm in a Single Copy Number Matrix
#'
#' This function iteratively applies the peelingOne function, thereby identifying multiple
#'
#' peaks across the genome in a single cohort.  Gains and losses should be analyzed separately.
#'
#' @param X A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
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
#' with mean copy number less than threshold will not be peeled.  Default = NULL.
#'
#' @param numIters The number of times peelingOne will be iterated.  Default = 5.
#'
#' @return A list containing two elements: X and interval.  X is an updated version of the input
#' copy number matrix in which the peak at k has been removed, and interval is genomic region
#' containing k.  By construction, interval cannot extend beyond the chromosome arm containing k.
#'
#' @details The \code{\link{peelingOne}} function applies the peeling algorithm described by 
#'   Walter et al. (Bioinformatics, 2011;27(5):678–685) at a given marker k. Because tumor genomes
#'   may contain multiple regions of copy number gain or loss, it important to apply the peeling
#'   process iteratively, thereby identifying multiple markers and surrounding genomic regions.
#' 
#' @examples lusc=pD[["X"]]
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
#'  peeledLusc=peelingOneIterate(X=lusc,posDT=posDT,gain=TRUE,nullDist=NULL,threshold=NULL,numIters=3)   
#'
#' @export

peelingOneIterate = function (X, posDT, gain = TRUE, nullDist = NULL, threshold = NULL, numIters = 5)
	{
	output = data.frame(chr = rep(NA, numIters), peelStart = rep(NA, numIters), peelStop = rep(NA, numIters),
		peakLoc = rep(NA, numIters), peakVal = rep(NA, numIters), peakPVal = rep(NA, numIters))
	gainLossInd = (-1)^(1 + gain)
	X1 = X * gainLossInd
	for (i in 1:numIters)
		{
		tempK = which.max(rowMeans(X1))
		tempPeakVal = gainLossInd * signif(rowMeans(X1)[tempK], 4)
		tempPeakPVal = ifelse(is.null(nullDist), NA, (1 + sum(nullDist > tempPeakVal))/length(nullDist) * (gainLossInd == 1) +
			(1 + sum(nullDist < tempPeakVal))/length(nullDist) * (gainLossInd == -1))
		tempPeakPval = min(tempPeakPVal, 1)
		tempPeel = peelingOne(X = X1, posDT = posDT, k = tempK, threshold = threshold)
		output[i,] = c(posDT[tempK, 1], tempPeel[[2]][1], tempPeel[[2]][2], posDT[tempK, 2], 
			tempPeakVal, tempPeakPVal)
		X1 = tempPeel[[1]]
		}
	return(output)
	}

