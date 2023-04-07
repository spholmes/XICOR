
#' Take fractionary part and make its binary expansion

#' Auxiliary function used in expanding real numbers
#' @param x is a number between 0 and 1
#' @return Binary expansion of length 31 of the decimal input
#' @note this implementation uses the built-in function intToBits


fracbinary<-function(x){
  ##Using integers in R requires the input
  ##to intToBits to be smaller than 2^31-1
  ##may change but see .Machine$integer.max
  return(rev(as.integer(intToBits((2^31)*x)))[-1])
}

