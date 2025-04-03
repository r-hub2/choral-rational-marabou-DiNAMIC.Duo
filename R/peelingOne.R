#' A Function to Apply the Peeling Algorithm in a Single Copy Number Matrix
#'
#' This function applies the peeling algorithm described in Walter et al. (Bioinformatics, 2011;27(5):678–685)
#'
#' to remove a peak from a copy number data set and define a genomic interval of interest around
#'
#' the peak.
#'
#' @param X A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
#'
#' @param posDT A data frame containing genomic position information for the genes in X.
#'
#' @param k The location (row of X) containing the peak that will be peeled.
#'
#' @param threshold A tuning parameter that controls the size of the peeled region.  Rows of X
#'
#' with mean copy number less than threshold will not be peeled.
#'
#' @return A list containing two elements: X and interval.  X is an updated version of the input
#' copy number matrix in which the peak at k has been removed, and interval is genomic region
#' containing k.  By construction, interval cannot extend beyond the chromosome arm containing k.
#'
#' @details Tumor genomes often contain multiple DNA copy number alterations, e.g., amplifications
#'   or deletions. The locus that harbors the most extreme alteration, k, as evidenced by the maximum
#'   or minimum column mean, provides a point estimate for the location of an underlying driver gene. 
#'   Also, loci near k may be affected by the same underlying genomic event. The peeling procedure is
#'   applied to "nullify" entries in X that contribute to the alteration at k, thus making it possible
#'   to identify altered regions elsewhere in the genome. This function is called by 
#'   \code{\link{peelingOneIterate}}.
#' 
#' @examples lusc=pD[["X"]]
#'
#'  posDT=pD[["posDT"]]
#'
#'  kLusc=which.max(rowMeans(lusc))
#'
#'  peeledLusc=peelingOne(X=lusc,posDT=posDT,k=kLusc,threshold=NULL)   
#'
#' @export

peelingOne = function (X, posDT, k, threshold = NULL) 
	{
	X = t(X)
    	n = dim(X)[1]
    	m = dim(X)[2]
    	rowMeans = rowMeans(X)
    	grandMean = mean(rowMeans)
    	chr = posDT[k, 1]
    	chrStart = min(which(posDT[,1] == chr))
    	chrEnd = max(which(posDT[,1] == chr))
    	chrM = chrEnd - chrStart + 1
    	chrK = k - chrStart + 1
    	chrData = X[,chrStart:chrEnd]
    	    	
	if (is.null(threshold))
		{
		colIndicator = as.numeric(colMeans(chrData) > grandMean)
		} else colIndicator = as.numeric(colMeans(chrData) > grandMean) * as.numeric(colMeans(chrData) > threshold)
  	colIndicator = dinamic::recodeBinary(binary.vec = colIndicator, k = chrK)

	exceed = rep(1, n) %o% colIndicator
    	whichBand = substring(posDT$cytoband[k], 1, 1)
    	bandVec = substring(posDT$cytoband[chrStart:chrEnd], 1, 1)
    	bandIndVec = as.numeric(bandVec == whichBand)
    	bandMatrix = rep(1, n) %o% bandIndVec
    	exceed = exceed * bandMatrix
    	peakInterval = as.numeric(colSums(exceed) > 0)
    	chrLeftSide = min(which(peakInterval == 1))
    	chrRightSide = max(which(peakInterval == 1))
    	leftSide = chrLeftSide + chrStart - 1
    	rightSide = chrRightSide + chrStart - 1
    	columnMean = mean(chrData[,chrK])
    	exceedRunning = matrix(0, n, chrM)
    	exceedRunning[,chrK] = as.numeric(chrData[,chrK] > rowMeans)
    	if (chrK > 1)
		{
        	for (j in (chrK - 1):1)
			{
            	exceedRunning[,j] = exceedRunning[,(j + 1)] * 
                		exceed[, j] * as.numeric(chrData[,j] > rowMeans)
        		}
    		}
    	if (chrK < chrM) 
		{
        	for (j in (chrK + 1):chrM) 
			{
            	exceedRunning[,j] = exceedRunning[,(j - 1)] * 
                		exceed[,j] * as.numeric(chrData[, j] > rowMeans)
        		}
    		}
    	notImputed = chrData * (1 - exceedRunning)
    	toBeImputed = chrData * exceedRunning
    	numTerm1 = n * grandMean
    	numTerm2 = sum(notImputed[,chrK])
    	numTerm3 = sum(exceedRunning[,chrK] * rowMeans)
    	denTerm = sum(toBeImputed[,chrK]) - numTerm3
    	tau = (numTerm1 - numTerm2 - numTerm3)/denTerm
    	imputed = tau * (toBeImputed - (matrix(rowMeans, n, chrM) * 
        exceedRunning)) + (matrix(rowMeans, n, chrM) * exceedRunning)
    	finalMatrix = imputed + notImputed
    	X[,chrStart:chrEnd] = finalMatrix
    	outputList = list(X = t(X), interval = c(posDT[leftSide, 2], posDT[rightSide, 3]))
    	return(outputList)
	}

