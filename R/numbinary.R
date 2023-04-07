#' Computes the binary expansion of a number
#' 
#' If the argument x is a real number the decimal
#' portion is dropped. 

#' @param x is a real or integer number
#' @return the output is a binary vector of length 31



numbinary<-function(x){
  ## Uses binary expansion  intToBits  
  ## Check that x is within the maximum integer range
  if ( x> .Machine$integer.max)  { 
    stop("overflow on the input, too large for an integer")}
  x<-as.integer(x)
  return(as.integer(intToBits(x))[-32])
}
