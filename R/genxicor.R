#' Compute the generalized cross rank increment correlation coefficient gxi.
#'
#' This function computes the generalized xi coefficient between two matrices
#' xmat and ymat.
#' There is a limitation on the size of the matrices, for the time
#' being, xmat and ymat can only have 31 columns.
#' If they are wider than 31, there is the option of using a
#' dimension reduction technique to bring the number of columns down
#' to 31, the first 31 components are then used.
#' The function encodes the data using a binary expansion and
#' then calls xicor on the vectors, so some of the arguments
#' relevant for xicor can be specified, such as pvalue.
#'
#' @aliases genxicor
#' @param xmat Matrix of numeric values in the first argument.
#' @param ymat Matrix of numeric values in the second argument.
#' @param pvalue Whether or not to return the p-value of rejecting
#' independence, if TRUE the function also returns the standard deviation of
#' xi.
#' 
#' @param method If method = "asymptotic" the function returns P-values
#' computed by the asymptotic theory. If method = "permutation", a permutation
#' test with nperm permutations is employed to estimate the P-value. Usually,
#' there is no need for the permutation test. The asymptotic theory is good
#' enough.
#' @param nperm In the case of a permutation test, \code{nperm} is the number
#' of permutations to do.
#' @return In the case pvalue=FALSE, function returns the value of the genxi
#' coefficient.
#' In the case pvalue=TRUE is chosen, the function returns a list:
#' \describe{\item{xi}{The
#' value of the xi coefficient.}
#' \item{sd}{The standard deviation.}
#' \item{pval}{The test p-value.}
#' }
#' @note This version does not use a seed as argument, if reproducibility is an issue, set a seed before calling the function.
#' @author Sourav Chatterjee, Susan Holmes
#' @export
#' @references Chatterjee, S. (2022) <arXiv:2211.04702>
#' @keywords ~methods ~htest
#' @examples 
#'



genxicor = function (xmat, ymat) {
  xmat <- as.matrix(xmat)
  ymat <- as.matrix(ymat)
  #check missingness in xmat and ymat
  if (any(is.na(xmat) ) || any(is.na(ymat)) )
  stop("One of the matrices xmat or ymat has missing values")
  n <- nrow(xmat)
  if (nrow(ymat) !=n)
    stop(" xmat and ymat need to have the same number of rows")
  x1 = rep(0,n)
  y1 = rep(0,n)
  for (i in 1:n) {
    x1[i] = as.numeric(borelmergeB(xmat[i,]))
    y1[i] = as.numeric(borelmergeB(ymat[i,]))
  }
  q = xicor(x1,y1, pvalue = TRUE)
  return(q)
}
