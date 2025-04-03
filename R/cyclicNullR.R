#' Create a cyclic shift-based null distribution for one or two copy number matrices
#'
#'
#' @param X a matrix or a data frame of copy number data. The rows and columns
#'   of X correspond to genes and subjects, respectively.
#'
#' @param Y a matrix or a data frame of copy number data. The rows and columns
#'   of X correspond to genes and subjects, respectively. It is assumed that
#'	 the rows of X and Y are indexed by the same set of genes that appear in
#'   genomic order. 
#'
#' @param numPerms the number of cyclic shifts used to create the null distribution.
#'   Default = 1e2.
#'
#' @param randomSeed a random seed.  Default = NULL.
#'
#' @return a matrix with two columns.  The first column, maxNull, is an empirical
#'   permutation-based null distribution for the maximum difference of row means 
#'   of X - row means of Y based on cyclic shift permutations of the columns of each 
#'   matrix; the second column, minNull, is an empirical distribution of the minimum 
#'   difference of the row means of X - the row means of Y based on the same permutations.
#'   If Y = NULL, the null distributions apply to the maximum and minimum row means of X.
#'
#' @details This function iteratively calls \code{\link{cyclicShiftColR}} to create an empirical
#'   permutation-based null distribution to assess the statistical significance of either
#'   (i) the maximum and minimum difference row means of X - row means of Y, or (ii) the 
#'   maximum and minimum row means of X, depending on whether two or one copy number matrices
#'   are being analyzed. The application of cyclic shift permutations to DNA copy number
#'   matrices was originally described by Walter et al. (Bioinformatics, 2011;27(5):678–685).
#'
#'
#' @examples
#' data(DiNAMIC.Duo)
#' output = cyclicNullR(X = pD[["X"]], Y = pD[["Y"]], numPerms = 25, randomSeed = NULL)
#'
#' @export

cyclicNullR = function(X, Y = NULL, numPerms = 1e2, randomSeed = NULL)
	{
	if (!is.null(Y))
		{
		if ((nrow(X) != nrow(Y)) || any(rownames(X) != rownames(Y)))
			{
			print("X and Y must have the same number of rows.  Also,")
			print("the rows of X and Y must be indexed by the same genes.")
			stop
			} else{  
				mX = ncol(X)
				mY = ncol(Y)
				maxNull = rep(NA, numPerms)
				minNull = rep(NA, numPerms)
				for (j in 1:numPerms)
					{
					Z = cbind(X, Y)
					if (!is.null(randomSeed))
						{
						Z = cyclicShiftColR(Z, randomSeed = j + randomSeed)
						} else Z = cyclicShiftColR(Z, randomSeed = NULL)
					maxNull[j] = max(rowMeans(Z[,1:mX]) - rowMeans(Z[,(mX + 1):(mX + mY)]))
					minNull[j] = min(rowMeans(Z[,1:mX]) - rowMeans(Z[,(mX + 1):(mX + mY)]))
					}
				}
			output = cbind(maxNull, minNull)
			colnames(output) = c("maxNullDist", "minNullDist")
		} else{
			maxNull = rep(NA, numPerms)
			minNull = rep(NA, numPerms)
			for (j in 1:numPerms)
				{
				Z = X
				if (!is.null(randomSeed))
					{
					Z = cyclicShiftColR(Z, randomSeed = j + randomSeed)
					} else Z = cyclicShiftColR(Z, randomSeed = NULL)
				maxNull[j] = max(rowMeans(Z))
				minNull[j] = min(rowMeans(Z))
				}
			output = cbind(maxNull, minNull)
			colnames(output) = c("maxNullDist", "minNullDist")
			}
	return(output)
	}

