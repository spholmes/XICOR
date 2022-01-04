#' Compute the FR coefficient  on two vectors based exactly on Gamma2.
#'
#' This function computes the unidimensional graph prediction coefficient
#' between two vectors xvec and yvec.
#'
#' @aliases  FRpredcor Gamma2
#' @param xvec Vector of numeric values in the first coordinate.
#' @param yvec Vector of numeric values in the second coordinate.
#' @param tiemethod Choice of treatment for ties, default is the "average"
#' @return In the case simple = TRUE, function returns the value of the
#'         FR standardized coefficient.
#' @note Auxiliary function with no checks for NA, etc.
#' @author Sourav Chatterjee, Susan Holmes
#' @seealso xicor FRpredcorhalf
#' @references Chatterjee, S. and Holmes, S (2020)
#' Practical observations and applications of the robust prediction
#' coefficient.
#' @keywords ~methods
#' @export
#' @examples
#' # Compute  the coefficient and compare to the xi coefficient
#' simulCompare <- function(n = 20, B = 1000)
#' {
#'  diffs<- rep(0,B)
#'  xvec <- 1:n
#'  for (i in 1:B)
#'  {
#'    yvec <- runif(n)
#'    diffs[i] <- FRpredcor(xvec, yvec) - xicor(xvec, yvec)
#'  }
#'  return(diffs)
#'  }
#'
#'  simulcompare1K <- simulCompare()
#'  summary(simulcompare1K)
#'
#'
#' @importFrom  stats complete.cases pnorm runif var

FRpredcor <- function(xvec, yvec, tiemethod= "average"){
  ### Two vectors, same length
  n <- length(xvec)
  ## Rearrange according to xvec
  PI <- rank(xvec, ties.method = tiemethod)
  ord <- order(PI)
  fr <- rank(yvec, ties.method = tiemethod)
  fr <- yvec[ord]
  R <- matrix(rep(fr,n), nrow = n, byrow = TRUE) - fr
  Rrank <- apply(abs(R), 1, rank, ties.method = tiemethod)
  FrRtotal <- 2 * sum(Rrank[row(Rrank)==(col(Rrank) - 1)])
return(FrRtotal)
}




