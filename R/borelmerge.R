#' Auxiliary function that takes avector and produces a single number through 
#' a Borel isomorphism using the wholebinary and backdec functions.

#' @param x is a vector of real numbers
#' @return produces a single real number by converting each element
# of x to binary, interlacing the digits, and reconverting it back to a single
# real number. (This is the Borel isomorphism.)



borelmerge = function(xvec) {
  
  n = length(xvec)
  rmatout = matrix(0, nrow = 2*n, ncol = 31)
  sgns = rep(0,n)
  
  for (i in 1:n) {
    w = wholebinary(xvec[i])
    ##Fill in the r matrix with two
    rmatout[2*i-1, ] = w$r[1,]
    rmatout[2*i, ] = w$r[2,]
    sgns[i] = w$s
  }
  if (ncol(rmatout)<32)
  {numberout<-backdec(weave(rmatout,sgns), 1)}
  else
  {stop("Need help here")}
  return(numberout)
}
