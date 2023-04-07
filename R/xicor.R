#' Compute the cross rank increment correlation coefficient xi.
#'
#' This function computes the xi coefficient between two vectors x and y,
#' possibly all coefficients for a matrix. If only one coefficient is computed
#' it can be used to test independence using a Monte Carlo permutation test or
#' through an asymptotic approximation test.
#'
#'
#' @aliases xi xicor
#' @param x Vector of numeric values in the first coordinate.
#' @param y Vector of numeric values in the second coordinate.
#' @param pvalue Whether or not to return the p-value of rejecting
#' independence, if TRUE the function also returns the standard deviation of
#' xi.
#' @param ties Do we need to handle ties? If ties=TRUE the algorithm assumes
#' that the data has ties and employs the more elaborated theory for
#' calculating s.d. and P-value. Otherwise, it uses the simpler theory. There
#' is no harm in putting ties = TRUE even if there are no ties.
#' @param method If method = "asymptotic" the function returns P-values
#' computed by the asymptotic theory. If method = "permutation", a permutation
#' test with nperm permutations is employed to estimate the P-value. Usually,
#' there is no need for the permutation test. The asymptotic theory is good
#' enough.
#' @param nperm In the case of a permutation test, \code{nperm} is the number
#' of permutations to do.
#' @param factor Whether to transform integers into factors, the default is to
#' leave them alone.
#' @return In the case pvalue=FALSE, function returns the value of the xi
#' coefficient, if the input is a matrix, a matrix of coefficients is returned.
#' In the case pvalue=TRUE is chosen, the function returns a list:
#' \describe{\item{xi}{The
#' value of the xi coefficient.}
#' \item{sd}{The standard deviation.}
#' \item{pval}{The test p-value.}
#' }
#' @note Dataset peas no longer available in psych, we are now using psychTools. 
#' @note This version does not use a seed as argument, if reproducibility is an issue, set a seed before calling the function.
#' @author Sourav Chatterjee, Susan Holmes
#' @seealso dcov
#' @import psychTools
#' @export
#' @references Chatterjee, S. (2020) <arXiv:1909.10140>.
#' @keywords ~methods ~htest
#' @examples
#'
#' ##---- Should be DIRECTLY executable !! ----
#' library("psychTools")
#' data(peas)
#' # Visualize       the peas data
#' library(ggplot2)
#' ggplot(peas,aes(parent,child)) +
#' geom_count() + scale_radius(range=c(0,5)) +
#'        xlim(c(13.5,24))+ylim(c(13.5,24))+       coord_fixed() +
#'        theme(legend.position="bottom")
#' # Compute one of the coefficients
#' xicor(peas$parent,peas$child,pvalue=TRUE)
#' xicor(peas$child,peas$parent)
#' # Compute all the coefficients
#' xicor(peas)
#'



xicor <- function(x, y = NULL, pvalue = FALSE, ties = TRUE, method = "asymptotic",
                                                                  nperm = 1000, factor = FALSE) {
  # x and y are the data vectors or a x is a matrix and y is null
  # to imitate the behavior of the standard cor function in R
	# The standard deviation of xi and the P-value for the 
  # test of independence is returned  if pvalue = TRUE. 
  # Otherwise, only the coefficient is returned.
	# If ties = TRUE, the algorithm assumes that the data has ties and employs the 
  # more elaborated theory for calculating s.d. and P-value. Otherwise, it uses the simpler theory. 
  # There is no harm in putting ties = TRUE even if there are no ties.
	# method = "asymptotic" returns P-values computed by the asymptotic theory.
  # If method = "permutation", a permutation test with nperm permutations is employed to 
  # estimate the P-value. Usually, there is no need for the permutation test. 
  # The asymptotic theory is good enough.
  # nperm is the number of permutations for the permutation test, if needed.
	# na.rm = TRUE results in the algorithm looking for NAs and removing them from x and y.
  # If it is known that the data has no NAs, one can set na.rm = FALSE and save a little time.
	# factor = TRUE results in the algorithm checking whether x and y are factor 
  # variables and converting them to integers if they are. 
  # If it is known that the variables are numeric, a little bit of time can be saved by 
  # setting factor = FALSE.
  #                                                             
  # Factor variables are converted to integers here:
            if (factor == TRUE) {
                            if (!is.numeric(x)) x <- as.numeric(factor(x))
                            if (!is.numeric(y)) y <- as.numeric(factor(y))
                                }
                          if (is.data.frame(y))
                                y <- as.matrix(y)
                          if (is.data.frame(x))
                                x <- as.matrix(x)
             if (!is.matrix(x) && is.null(y))
                          stop("supply both 'x' and 'y' or a matrix-like 'x'")
             if (!(is.numeric(x) || is.logical(x)))
                          stop("'x' must be numeric")
             stopifnot(is.atomic(x))
             if (!is.null(y)) {
                          if (!(is.numeric(y) || is.logical(y)))
                                                stop("'y' must be numeric")
                          stopifnot(is.atomic(y))
             }
             if (is.null(y)) {
                          ncy <- ncx <- ncol(x)
                          if (ncx == 0)
                                    stop("'x' is empty")
                          if (pvalue == TRUE)
                                    stop("testing is not available for matrices")
              r <- matrix(0, nrow = ncx, ncol = ncy)
              for (i in seq_len(ncx)) {
                        for (j in seq_len(i)) {
                               x2 <- x[, i]
                               y2 <- x[, j]
                               ok <- complete.cases(x2, y2)
                               x2 <- x2[ok]
                               y2 <- y2[ok]
                               if (any(ok))
                               { r[i, j] <- calculateXI(x2, y2, simple=TRUE)
                               ###it's not symmetric, we have to compute both
                                 r[j, i] <- calculateXI(y2, x2, simple=TRUE)
                               }
                               else NA
              }
              }
              rownames(r) <- colnames(x)
              colnames(r) <- colnames(x)
              return(r)
       }
       ##Two vectors case
       else
              if (ncol(as.matrix(x))==1 & ncol(as.matrix(y))==1){
              ok <- complete.cases(x, y)
              x <- x[ok]
              y <- y[ok]
              res <- calculateXI(x, y, simple = FALSE)
              xi <- res$xi
              CU <- res$CU
              n <- length(x)
        }
       else
              if (ncol(as.matrix(x))>1 & ncol(as.matrix(y))==1){
                        ok <- complete.cases(cbind(x, y))
                        x <- x[ok]
                        y <- y[ok]
                        res <- calculateXI(x, y, simple = FALSE)
                        xi <- res$xi
                        CU <- res$CU
                        n <- length(x)
              }
	# If P-value needs to be computed:
	if (pvalue) {                         
		# If there are no ties, return xi, the s.d., and theoretical P-value:
		if (ties == FALSE) return(list(xi = xi, sd = sqrt(2/(5*n)),
                                           pval = 1 - pnorm(sqrt(n)*xi/sqrt(2/5))))
	       if (!(method %in% c("asymptotic","permutation")))
	       stop("method for test can only be asymptotic or permutation")
	       # If there are ties, and the theoretical method is to be used
	       if (method == "asymptotic") {
		       fr <- res$fr
			# The following steps calculate the theoretical variance in the presence of ties:
			qfr <- sort(fr)
			ind <- c(1:n)
			ind2 <- 2*n - 2*ind + 1
			ai <- mean(ind2*qfr*qfr)/n
			ci <- mean(ind2*qfr)/n
			cq <- cumsum(qfr)
			m <- (cq + (n - ind)*qfr)/n
			b <- mean(m^2)
			v <- (ai - 2*b + ci^2)/(CU^2)
			# Return xi, standard deviation of xi, and P-value:
			return(list(xi = xi, sd = sqrt(v/n), pval = 1 - pnorm(sqrt(n)*xi/sqrt(v))))
		}
    #
		# If permutation test is to be used for calculating P-value:
		if (method == "permutation") {
			rp <- rep(0, nperm)
			for (i in 1:nperm) {
				x1 <- runif(n, 0, 1)
				rp[i] <- calculateXI(x1, y)
			}
			# Return P-value and sd based on permutation test:
			return(list(xi = xi, sd = sqrt(var(rp)), pval = mean(rp > xi)))
		}
	}
	# If only xi is desired, return value for xi:
	else return(xi)
}             
