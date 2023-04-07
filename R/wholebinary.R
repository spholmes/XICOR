#' Encodes a number as a two row binary matrix and its sign
#' 
#' Auxiliary function used for generating expansion of a number,
#' the binary expansion of length nc of the integer part is the first row 
#' and the binary expansion of length nc of the fractional part
#' is the second row of the matrix.
#' The sign as appended into the final list object which the function
#' returns.
#' 
#' @param x is a decimal number 
#' @param nc is the length of the binary expansion and defines the
#' number of columns of the output matrix
#' @return This function generates a list with a binary matrix
#' rmat
#' with two rows and the sign sgn in a separate entry of the list.
wholebinary = function(x,nc=31) {
  x1 = abs(x)
  s = sign(x)
  if (s == 0) s = 1
  y = floor(x1)
  z = x1 - y
  rmatrix = matrix(nrow = 2, ncol = nc)
  rmatrix[1,] = numbinary(y)
  rmatrix[2,] = fracbinary(z)
  return(list(rmat = rmatrix, sgn = s))
}
