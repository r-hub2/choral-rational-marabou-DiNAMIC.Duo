#' Perform the cyclic shift procedure on the columns of a matrix
#'
#'
#' @param X a matrix or a data frame of copy number data. The rows and columns
#'   of X correspond to genes and subjects, respectively.
#'
#' @param randomSeed a random seed.  Default = NULL.
#'
#' @return a matrix Z whose dimensions are the same as X.  Each column of Z
#'   is obtained by perform a cyclic shift of the corresponding column of X.
#'
#' @details Many algorithms for identifying recurrent DNA copy number alterations,
#'   e.g., amplifications and deletions, assess statistical significance using
#'   permutation-based null distributions.  Like many genomic data types, DNA copy 
#'   number data has an underlying spatial correlation structure.  Randomly permuting
#'   the values within a given subject ignores this structure, and this can impact
#'   type I error rates.  In contrast, cyclic permutation largely preserves the
#'   local correlation structure.
#'
#'
#' @examples
#' test = matrix(c(1:50), 10, 5)
#' cyclicShiftColR(test, randomSeed = NULL)
#'
#' @export

cyclicShiftColR = function(X, randomSeed = NULL)
	{
	n = nrow(X)
	m = ncol(X)
	if (!is.null(randomSeed))
		{
		set.seed(randomSeed)
		}
	Z = X
	shiftVec = sample(1:n, m, replace = TRUE)
	for (i in 1:m)
		{
		k = shiftVec[i]
		Z[,i] = ((k == 1) * Z[,i]) + ((k != 1) * Z[c(c(k:n), c(1:(k - 1))), i][c(1:n)])
		}
	return(Z)	
	}

