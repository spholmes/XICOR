#' Inverse function to wholebinary returns the number from its expansion
#' 
#' @param rmat is a matrix of two rows, the first row
#' of the matrix is the expansion of the integer part
#' the second row is the binary expansion of the fractional part.
#' @param sgn is the sign
#' @note It may be necessary to make a new version of this
#' using special functions for large integers.

backdec = function(rmat,sgn) {
  # Need a backdecx function that supports large integers  
  p = ncol(rmat)
  q1 = 2^c(0:(p-1))
  q2 = 2^(-c(1:p))
  return(sgn*(sum(rmat[1,]*q1) + sum(rmat[2,]*q2)))
}

