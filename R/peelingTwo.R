#' A Function to Apply the Peeling Algorithm in a Two Copy Number Matrices
#'
#' This function applies a modified version of the peeling algorithm originally described in Walter et al.,
#'
#' (PMID  21183584) to remove a peak from the copy number differences and define a genomic interval of interest
#'
#' around the peak.
#'
#' @param X A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
#'
#' @param Y A matrix of normalized gene-level copy number data (rows = genes, columns = subjects).
#'
#' @param posDT A data frame containing genomic position information for the genes in X.
#'
#' @param k The location (row of X and Y) containing the peak that will be peeled.
#'
#' @param threshold A tuning parameter that controls the size of the peeled region.  Rows in which 
#'
#' rowMeans(X) - rowMeans(Y) are less than threshold will not be peeled.
#'
#' @return A list containing three elements: X, Y, and interval.  X and Y are updated versions of the
#' input copy number matrices X and Y in which the peak at k has been removed, and interval is genomic region
#' containing k.  By construction, interval cannot extend beyond the chromosome arm containing k.
#'
#' @details When tumor genomes from two cohorts are compared, there may be multiple regions that harbor
#'   copy number differences.  For example, gains or losses may be present in only one of the two cohorts,
#'   and this could give rise to copy number differences.  Alternatively, the same region of the genome
#'   may exhibit gain or loss in both cohorts.  If the magnitudes of the common gain or loss are distinct,
#'   then this also gives rise to copy number differences.  The locus that harbors the most extreme difference,
#'   k, provides a point estimate for the underlying driver gene that gives rise to the difference.  Loci
#'   near k may also be affected by the underlying difference in copy number.   The peeling procedure for
#'   two cohorts is applied to "nullify" entries of both X and Y that contribute to the alteration at k,
#'   thus making it possible to identify other regions of the genome that harbor copy number differences.
#'   This function is called by \code{\link{peelingTwoIterate}}.
#' 
#' @examples luad=pD[["X"]]
#'
#'  lusc=pD[["Y"]]
#'
#'  posDT=pD[["posDT"]]
#'
#'  kDiff=which.max(rowMeans(luad)-rowMeans(lusc))
#'
#'  peeledDiff=peelingTwo(X=luad,Y=lusc,posDT=posDT,k=kDiff,threshold=NULL)
#'
#' @export

peelingTwo = function(X, Y, posDT, k, threshold = NULL)
	{
	X = t(X)
	Y = t(Y)
  	nX = dim(X)[1]
	m = dim(X)[2] # x: nx by mx (nx = number of cases in X, m = number of markers)
  	nY = dim(Y)[1]
	mY = dim(Y)[2] # y: ny by my
  	if (!m==mY)
		{
    		stop("The number of markers in X and Y should be the same")
  		}
  	rowMeansX = rowMeans(X, na.rm=T) 
  	rowMeansY = rowMeans(Y, na.rm=T) 
  	grandMeanX = mean(rowMeansX) 
  	grandMeanY = mean(rowMeansY)

  	#This section of code creates exceed, a binary n x chr.m matrix, where chr.m is the number of markers
  	#in the chromosome containing k, the marker where the peeling procedure begins.  An entry of exceed is 1
  	#if the corresponding entry of x belongs to the peak interval
	chr = posDT[k, 1]
    	chrStart = min(which(posDT[,1] == chr))
    	chrEnd = max(which(posDT[,1] == chr))
    	chrM = chrEnd - chrStart + 1
    	chrK = k - chrStart + 1
    	chrDataX = X[,chrStart:chrEnd]
	chrDataY = Y[,chrStart:chrEnd]

  	colMeansX = colMeans(chrDataX)
  	colMeansY = colMeans(chrDataY)
  	
  	if (is.null(threshold))
		{
		colIndicator = as.numeric((colMeansX - colMeansY) > (grandMeanX - grandMeanY))
		} 	else colIndicator = as.numeric((colMeansX - colMeansY) > (grandMeanX - grandMeanY)) * 
				as.numeric((colMeansX - colMeansY) > threshold)

  	colIndicator = dinamic::recodeBinary(binary.vec = colIndicator, k = chrK)

	whichBand = substring(posDT$cytoband[k], 1, 1)
    	bandVec = substring(posDT$cytoband[chrStart:chrEnd], 1, 1)
    	bandIndVec = as.numeric(bandVec == whichBand)

  	exceedX = rep(1, nX) %o% colIndicator
  	bandMatrixX = rep(1, nX) %o% bandIndVec
  	exceedX = exceedX * bandMatrixX
  
  	exceedY = rep(1, nY) %o% colIndicator
  	bandMatrixY = rep(1, nY) %o% bandIndVec
  	exceedY = exceedY * bandMatrixY
  
  	#This section of code find the markers that make up the peak interval
  	peakInterval = as.numeric(colSums(exceedX) > 0)
  	chrLeftSide = min(which(peakInterval == 1))
  	chrRightSide = max(which(peakInterval == 1))
  	leftSide = chrLeftSide + chrStart - 1
  	rightSide = chrRightSide + chrStart - 1

  	#These will be used below
  	columnMeanX = mean(chrDataX[,chrK])
  	columnMeanY = mean(chrDataY[,chrK])
  
  	#This section of code creates exceed.running, a binary n x chr.m matrix.  An entry of exceed.running
  	#is 1 if the corresponding entry of x will be peeled, otherwise it is zero.
  	exceedRunningX = matrix(0, nX, chrM)
  	exceedRunningX[,chrK] = as.numeric(chrDataX[,chrK] > rowMeansX)
  	if (chrK > 1)
  		{
    		for (j in (chrK - 1):1)
    			{
      		exceedRunningX[,j] = exceedRunningX[,(j + 1)] * exceedX[,j] * as.numeric(chrDataX[,j] > rowMeansX)
    			}
  		}
  	if (chrK < chrM)
  		{
    		for (j in (chrK + 1):chrM)
    			{
      		exceedRunningX[,j] = exceedRunningX[,(j - 1)] * exceedX[,j] * as.numeric(chrDataX[,j] > rowMeansX)
    			}
  		}

  	exceedRunningY = matrix(0, nY, chrM)
  	exceedRunningY[,chrK] = as.numeric(chrDataY[,chrK] < rowMeansY)
  	if (chrK > 1)
  		{
    		for (j in (chrK - 1):1)
    			{
      		exceedRunningY[,j] = exceedRunningY[,(j + 1)] * exceedY[,j] * as.numeric(chrDataY[,j] < rowMeansY)
    			}
  		}
  	if (chrK < chrM)
  		{
    		for (j in (chrK + 1):chrM)
    			{
      		exceedRunningY[,j] = exceedRunningY[,(j - 1)] * exceedY[,j] * as.numeric(chrDataY[,j] < rowMeansY)
    			}
  		}
  
  	#This section of code implements a multiplicative-based shrinkage approach at column
  	ixSub = which(chrDataX[,chrK] > rowMeansX)
  	iySub = which(chrDataY[,chrK] < rowMeansY)
  
  	A = sum(chrDataX[ixSub,chrK] - rowMeansX[ixSub])/nX
  	B = sum(chrDataY[iySub,chrK] - rowMeansY[iySub])/nY
  	C = (sum(rowMeansX[ixSub]) + sum(chrDataX[(!1:nX %in% ixSub),chrK]))/nX
  	D = (sum(rowMeansY[iySub]) + sum(chrDataY[(!1:nY %in% iySub),chrK]))/nY
  
  	tau = (grandMeanX - grandMeanY - C + D)/(A - B)
  	finalMatrixX = (tau * (chrDataX - (rowMeansX %o% rep(1, ncol(chrDataX)))) + (rowMeansX %o% rep(1, ncol(chrDataX)))) * exceedRunningX + 
		chrDataX * (1 - exceedRunningX)
  	finalMatrixY = (tau * (chrDataY - (rowMeansY %o% rep(1, ncol(chrDataY)))) + (rowMeansY %o% rep(1, ncol(chrDataY)))) * exceedRunningY + 
		chrDataY * (1 - exceedRunningY)
  
  	#Create the output
  	X[,chrStart:chrEnd] = finalMatrixX
  	Y[,chrStart:chrEnd] = finalMatrixY
  	outputList = list(X = t(X), Y = t(Y), interval = c(posDT[leftSide, 2], posDT[rightSide, 3]))
  
  	#Return the output
  	return(outputList)
	}

