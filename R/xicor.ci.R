#' Computes confidence intervals for rank correlation coefficient xi.
#'
#' This function implements different methods for computing a confidence
#' interval for xicor. As the common n-out-of-n bootstrap has been shown
#' to fail both for estimating the variance [1] and for estimating the
#' distribution of xicor [3], bootstrap methods utilize the m-out-of-n
#' bootstrap. The parameter `methods` can be any combination of the
#' following strings, or "all" for all methods:
#' \describe{\item{"boot.norm"}{An interval based on asymptotic normality
#' and root-n consistency of xicor. The standard deviation is estimated
#' with an n-choose-m bootstrap. Suggested by Dette & Kroll [2].}
#' \item{"boot.nonpar"}{A non-parametric n-choose-m bootstrap interval
#' based on root-n consistency of xicor. The bootstrap is done with a
#' normalized version of xicor that has lower bias. Suggested by
#' Dalitz, Arning, & Goebbels [3].}
#' \item{"sigma.lh"}{An interval based on asymptotic normality for
#' continuous x and y. The standard deviation is estimated with the
#' second estimator (overlined sigma) by Lin & Han [4]. As this estimator
#' can yield negative values for the variance, especially for small sample
#' sizes, computation of a confidence interval with this method can fail.}
#' }
#'
#' @param x Vector of numeric values in the first coordinate.
#' @param y Vector of numeric values in the second coordinate.
#' @param conf The confidence level.
#' @param method The method for constructing the interval; can be any combination of
#'        "boot.norm", "boot.nonpar", "sigma.lh", or "all" for all methods.(see Description).
#' @param R The number of bootstrap repetitions.
#' @return A data frame with columns:
#'   \item{method}{The method used for the confidence interval}
#'   \item{lower}{The lower bound of the confidence interval}
#'   \item{upper}{The upper bound of the confidence interval}
#' 
#'   The confidence level is stored as an attribute "conf".
#' @author Christoph Dalitz, Juliane Arning, Felix Lögler
#' @seealso xicor
#' @export
#' 
#' @references 
#' \describe{\item{[1]}{Lin, Z. and Han, F. (2024):
#' _On the failure of the bootstrap for Chatterjee’s rank correlation_
#' Biometrika, asae004, 02,
#' <doi:10.1093/biomet/asae004>.}
#' \item{[2]}{Dette, H. and Kroll, M. (2024):
#' _A simple bootstrap for Chatterjee’s rank correlation_
#' Biometrika, asae045,
#' <doi:10.1093/biomet/asae045>.}
#' \item{[3]}{Dalitz, C., Arning, J. and Goebbels, S. (2024)
#' _A Simple Bias Reduction for Chatterjee’s Correlation_
#' Journal of Statistical Theory and Practice 18, 51,
#' <doi:10.1007/s42519-024-00399-y>.}
#' \item{[4]}{Lin, H. and Han, F. (2022):
#' _Limit theorems of Chatterjee’s rank correlation_
#' <arxiv:2204.08031>.}
#' }
#' 
#' 
#' @examples
#'
#' x <- runif(100, min=-1, max=1)
#' y <- x^2 + rnorm(100, sd=0.2)
#' ci <- xicor.ci(x, y)
#' print(ci)
#'
#' # Compute CI using only bootstrap normal method
#' ci_boot_norm <- xicor.ci(x, y, method = "boot.norm")
#' print(ci_boot_norm)
#' 

xicor.ci <- function(x, y, conf = 0.95, method = "all", R = 2000) {

  allmethods <- c("boot.norm", "boot.nonpar", "sigma.lh")
  if (method == "all") method <- allmethods

  # some plausi checks
  stopifnot(length(x) == length(y))
  stopifnot(0 < conf & conf < 1)
  stopifnot(all(method %in% allmethods))

  # allocate result data.frame
  res <- data.frame(method = method, lower = 0, upper = 0)
  attr(res, "conf") <- conf

  # apply reqeusted methods
  if ("boot.nonpar" %in% method) {
    res[res$method == "boot.nonpar",] <- xicor.ci.boot.nonpar(x, y, conf, R)
  }
  if ("boot.norm" %in% method) {
    res[res$method == "boot.norm",] <- xicor.ci.boot.norm(x, y, conf, R)
  }
  if ("sigma.lh" %in% method) {
    res[res$method == "sigma.lh",] <- xicor.ci.sigma.lh(x, y, conf)
  }

  return(res)
}


#
# actual implementation of the methods
#

#' Bootstrap normality method for xicor confidence interval
#' @keywords internal
#' @noRd
#' 
xicor.ci.boot.norm <- function(x, y, conf = 0.95, R = 2000) {

  # wrapper around xicor() for bootstrap
  xiboot <- function(indices, data) {
    xicor(data[indices, "x"], data[indices, "y"])
  }

  # use a matrix instead of a data frame to speed up apply()
  data <- as.matrix(cbind(x = x, y = y))
  n <- nrow(data)

  # point estimate
  xi <- xiboot(1:n, data)

  # bootstrap distribution with m chosen according to Dette & Kroll (2023)
  m <- floor(n^(1 / 2))
  indices <- replicate(R, sample(1:n, size = m, replace = F))
  xi.star <- apply(indices, MAR = 2, xiboot, data = data)

  # normal distribution confidence interval with correction factor tau
  tau <- sqrt
  s <- tau(m) * sd(xi.star, na.rm = T)
  z <- qnorm((1 + conf) / 2)
  return(data.frame(method = "boot.norm",
                    lower = xi - z * s / tau(n),
                    upper = xi + z * s / tau(n)))
}


#' Nonparametric method for xicor confidence interval
#' @keywords internal
#' @noRd
#' 
xicor.ci.boot.nonpar <- function(x, y, conf = 0.95, R = 2000) {

  # wrapper around xicor() for bootstrap
  # xicor is normalized for better coverage probability
  xiboot <- function(indices, data) {
    xicor(data[indices, "x"], data[indices, "y"]) /
      xicor(data[indices, "y"], data[indices, "y"])
  }

  # use a matrix instead of a data frame to speed up apply()
  data <- as.matrix(cbind(x = x, y = y))
  n <- nrow(data)

  # point estimate
  xi <- xiboot(1:n, data)

  # bootstrap distribution with m chosen according to Dalitz et al. (2023)
  m <- round(2 * sqrt(n))
  indices <- replicate(R, sample(1:n, size = m, replace = F))
  xi.star <- apply(indices, MAR = 2, xiboot, data = data)

  # confidence interval according to Politis & Romano (1994)
  tau <- sqrt
  xq <- tau(m) * (quantile(xi.star, c((1 + conf) / 2, (1 - conf) / 2), na.rm = TRUE) - xi)
  return(data.frame(method = "boot.nonpar",
                    lower = xi - xq[1] / tau(n),
                    upper = xi - xq[2] / tau(n)))
}


#' Implementation of method "sigma.bar"
#' @keywords internal
#' @noRd
#' 
xicor.ci.sigma.lh <- function(x, y, conf = 0.95) {
  if (length(x) != length(unique(x)) | length(y) != length(unique(y)))
    warning("method 'sigma.lh' is meant for continuous variables")

  # sort (x,y) after x for simple indexing (right) NN of X
  order_index <- order(x)
  x <- x[order_index]
  y <- y[order_index]

  R <- rank(y, ties.method = "max") # R denote the rank of Y
  n <- length(x)

  N <- matrix(nrow = 3, ncol = n) # to simplify Nbar is refered to as N
  for (i in 1:n) {
    # set the NN or i if i is under the k largest
    N[1, i] <- ifelse(i >= n, i, i + 1)
    N[2, i] <- ifelse(i >= n - 1, i, i + 2)
    N[3, i] <- ifelse(i >= n - 2, i, i + 3)
  }

  # Precalculate min(R[i], R[N(i)]) and similar
  mins.ri.rni <- pmin(R, R[N[1,]])
  mins.ri.rn2i <- pmin(R, R[N[2,]])
  mins.rn2i.rn3i <- pmin(R[N[2,]], R[N[3,]])
  mins.rni.rn2i <- pmin(R[N[1,]], R[N[2,]])

  t1 <- 1 / n^3 * sum(mins.ri.rni^2)
  t2 <- 2 / n^3 * sum(mins.ri.rni * mins.ri.rn2i)
  t3 <- -2 / n^3 * sum(mins.ri.rni * mins.rn2i.rn3i)

  # R[i] <= min(R[j], R[N[1, j]]) for each i,j j!=i
  t45.outer <- outer(R, mins.ri.rni, "<=")

  t4.help <- t45.outer * rep(mins.ri.rni, times = n)
  diag(t4.help) <- 0 #
  t4 <- 4 / (n^2 * (n - 1)) * sum(t4.help)

  t5.help <- t45.outer * rep(mins.rni.rn2i, times = n)
  diag(t5.help) <- 0 # exclude i=j
  t5 <- -2 / (n^2 * (n - 1)) * sum(t5.help)

  t6.help <- outer(mins.ri.rni, mins.ri.rni, pmin)
  diag(t6.help) <- 0
  t6 <- 1 / (n^2 * (n - 1)) * sum(t6.help)

  t7 <- -4 * (1 / n^2 * sum(mins.ri.rni))^2
  sum <- t1 + t2 + t3 + t4 + t5 + t6 + t7

  sigma <- ifelse(sum >= 0, sqrt(36*sum), NA)
  xi <- xicor(x, y)
  z <- qnorm((1 + conf) / 2)
  return(data.frame(method = "sigma.lh",
                    lower = xi - z * sigma / sqrt(n),
                    upper = xi + z * sigma / sqrt(n)))
}
