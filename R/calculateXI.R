#' Compute the cross rank coefficient xi on two vectors.
#'
#' This function computes the xi coefficient between two vectors x and y.
#'
#' @aliases  xicorcoefficient
#' @param xvec Vector of numeric values in the first coordinate.
#' @param yvec Vector of numeric values in the second coordinate.
#' @param simple Whether auxiliary information is kept to pass on.
#' @return In the case simple = TRUE, function returns the value of the
#' xi
#' coefficient,
#' If simple = FALSE is chosen, the function returns a list:
#' \describe{\item{xi}{The xi coefficient}
#' \item{fr}{rearranged rank of yvec}
#' \item{CU}{\code{mean(gr*(1-gr))}}
#' }
#' @note Auxiliary function with no checks for NA, etc.
#' @author Sourav Chatterjee, Susan Holmes
#' @seealso xicor
#' @references Chatterjee, S. (2020) A New Coefficient Of Correlation,
#' <arXiv:1909.10140>.
#' @keywords ~methods
#' @export
#' @examples
#' # Compute one of the coefficients
#' library("psychTools")
#' data(peas)
#' calculateXI(peas$parent,peas$child)
#' calculateXI(peas$child,peas$parent)


calculateXI <- function(xvec, yvec, simple=TRUE) {
## The following function computes the new correlation coefficient.
## Main simple correlation calculation in the case of two vectors xvec and yvec, no missing.
## This will eventually benefit from being written in C for speed.
## This version does not have a seed that can be fixed.
# n is the sample size.
  n <- length(xvec)							
# PI is the rank vector for x, with ties broken at random
  PI <- rank(xvec, ties.method = "random")	
# fr[i] is number of j s.t. y[j] <= y[i], divided by n.  
  fr <- rank(yvec, ties.method = "max")/n		
# gr[i] is number of j s.t. y[j] >= y[i], divided by n.  
  gr <- rank((- yvec), ties.method = "max")/n	
# order of the x's, ties broken at random.
  ord <- order(PI)					
# Rearrange fr according to ord.  
  fr <- fr[ord]								
# xi is calculated in the next three lines:
  A1 <- sum(abs(fr [1:(n - 1)] - fr [2:n]))/ (2*n)
  CU <- mean(gr* (1 - gr))
  xi <- 1 - A1/CU
  if (simple == TRUE)
    return(xi)
  else
    return(list(xi = xi,fr = fr,CU = CU))
}
