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

#' @return Function returns the value of the genxi coefficient.
#' Since by default the option pvalue=TRUE is chosen, the function 
#' returns a list:
#' \describe{\item{xi}{The
#' value of the xi coefficient.}
#' \item{sd}{The standard deviation.}
#' \item{pval}{The test p-value.}
#' }
#' @note This version does not use a seed as argument, if reproducibility is an issue, set a seed before calling the function.
#' @note The p-value of rejecting independence is set to TRUE.
#' @author Sourav Chatterjee, Susan Holmes
#' @export
#' @references Chatterjee, S. (2022) <arXiv:2211.04702>
#' @keywords ~methods ~htest
#' @examples 
#' 
#' example_joint_calc = function(n,x=runif(n),y=runif(n),ep=runif(n)) {
#' u = (x + y + ep) %% 1
#' v = ((x + y)/2 + ep) %% 1
#' w = (4*x/3 + 2*y/3 + ep) %% 1
#' z = (2*x/3 + y/3 + ep) %% 1
#' q = cbind(u,v,w,z)
#' p = cbind(x,y)
#' c1 = genxicor(u, p)
#' c2 = genxicor(v, p)
#' c3 = genxicor(w, p)
#' c4 = genxicor(z, p)
#' c5 = genxicor(q, p)
#' return(list(marg1 = c1$xi, marg2 = c2$xi, marg3 = c3$xi, 
#' marg4 = c4$xi, joint = c5$xi, p1 = c1$pval, p2 = c2$pval, p3 = c3$pval,
#' p4 = c4$pval, p5 = c5$pval))
#' }
#' result1 <- example_joint_calc(n=10)
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
    x1[i] = as.numeric(borelmerge(xmat[i,]))
    y1[i] = as.numeric(borelmerge(ymat[i,]))
  }
  q = xicor(x1,y1, pvalue = TRUE)
  return(q)
}
