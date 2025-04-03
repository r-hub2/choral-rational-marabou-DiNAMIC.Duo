#' Find DiNAMIC's null distribution
#'
#'
#' @param x An n by m numeric matrix containing DNA copy number data from n subjects at m markers.
#'
#' @param num.perms A positive integer that represents the number of cyclic shifts used to create 
#'   the empirical distribution.
#'
#' @param random.seed An optional random seed (default = NULL).
#' 
#' @return A numerical vector of length num.perms
#'
#' @details The cyclic shift procedure is detailed in Bioinformatics (2011) 27(5) 678 - 685.  Briefly, 
#'   cyclic shift is a permutation procedure for DNA copy number data that largely preserves the underlying 
#'   correlation of the markers.  This function uses \code{num.perms} cyclic shifts of the copy number matrix 
#'   \code{x} to create an approximate null distribution for \code{max(colSums(x))} or \code{min(colSums(x))}.  
#'   The statistical significance of the observed value of \code{max(colSums(x))} or \code{min(colSums(x))} 
#'   is assessed by the functions \code{\link{quickLook}} and \code{\link{detailedLook}}.
#'
#'
#' @examples
#' random.seed = 12345
#' set.seed(random.seed)
#' x = matrix(rnorm(50), 5, 10)
#' num.perms = 10
#' findNull(x, num.perms, random.seed)
#'
#' @export

findNull = function(x, num.perms, random.seed = NULL) 
	{
    if (length(random.seed) > 0) {
        set.seed(random.seed)
		}
    n = dim(x)[1]
    m = dim(x)[2]
    shift.vec = rep(23, num.perms)
    for (j in (1:num.perms)) {
        cutpoints = sample(c(1:m), size = n, replace = TRUE)
        perm.x = matrix(0, n, m)
        for (i in (1:n)) {
            index = cutpoints[i]
            if (index == 1) {
                second = c()
            }
            if (index > 1) {
                second = x[i, 1:(index - 1)]
            }
            first = x[i, index:m]
            perm.x[i, ] = c(first, second)
			}
        perm.col.sums = colSums(perm.x)
        shift.vec[j] = max(colSums(perm.x))
    }
    return(shift.vec)
	}
