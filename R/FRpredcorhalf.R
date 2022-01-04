#' Compute the FR half coefficient  on two vectors based on half Gamma 2.
#'
#' This function computes the unidimensional ranked
#' half graph prediction coefficient
#' between two vectors xvec and yvec.
#'
#' @aliases  FRpredcorhalf
#' @param xvec Vector of numeric values in the first coordinate.
#' @param yvec Vector of numeric values in the second coordinate.
#' @param tiemethod Choice of treatment for ties, default is the "average"
#' @return In the case simple = TRUE, function returns the value of the
#'         FR standardized coefficient.
#' @note Auxiliary function with no checks for NA, etc.
#' @author Sourav Chatterjee, Susan Holmes
#' @seealso xicor FRpredcor
#' @references Chatterjee, S. and Holmes, S (2020)
#' Practical observations and applications of the robust prediction
#' coefficient.
#' @keywords ~methods
#' @export
#' @examples
#' # Compute  the coefficient and compare to the xi coefficient
#' simulCompare <- function(n = 20, B = 1000)
#' {
#'  diffsim <- rep(0,B)
#'  xvec <- 1:n
#'  for (i in 1:B)
#'  {
#'    yvec <- sample(n,n)
#'    diffsim[i] <- FRpredcorhalf(xvec,yvec)-xicor(xvec,yvec)
#'  }
#'  return(diffsim)
#'  }
#'
#'  compare1K <- simulCompare()
#'  summary(compare1K)
#'
#'
#' @importFrom  stats complete.cases pnorm runif var

FRpredcorhalf <-function(xvec, yvec, tiemethod= "average"){
  ### Two vectors, same length
  n <- length(xvec)
  ## Rearrange according to xvec
  PI <- rank(xvec, ties.method = tiemethod)
  ord <- order(PI)
  fr <- rank(yvec, ties.method = tiemethod)
  fr <- fr[ord]
  R <- matrix(rep(fr,n), nrow =n, byrow=TRUE) - fr
  FrR <- sum(diag(abs(R[1:(n - 1),2:n])))
  Frcoeff <- 1 - (3*FrR) / (n^2 - 1)
return(Frcoeff)
}







