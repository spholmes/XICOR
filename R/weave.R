#' Take a matrix of two numbers given in their binary expansion
#' one in each of the two rows
#' and return the interleaving of the two numbers
#' 
#' @param rmat a matrix with two times m rows corresponding to the 
#' the expansions of the m numbers to be interleaved.
#' @param sgn 
#' 
weave = function(rmat, sgn) {
  ##weave the two numbers and keep an extra
  ##rmat is a matrix coming from 
  ##sign element??  
  n = nrow(rmat)
  p = ncol(rmat)
  m = n/2
  r1 = matrix(0, nrow = 2, ncol = m*p + m + 1) 
  w = m*c(0:(p-1))
  for(i in 1:m) {
    r1[1,w+i] = rmat[2*i-1,]
    r1[2,w+i] = rmat[2*i,]
  } 
  sgn = (sgn + 1)/2
  l = which(r1[1,] == 1)
  if (sum(l) != 0) l1 = max(l)
  else l1 = 0
  r1[1,(l1+1):(l1+m)] = sgn
  r1[1, l1+m+1] = 1
  return(r1)
}