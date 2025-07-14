#' mvlk
#'
#' @param C A linear kernel in the form of a symmetric numeric matrix. Row-
#' and column-order must match the order of the levels of the environment
#' factor in the dataframe that will be passed to ASReml-R.
#'
#' @return A function compatible with asreml::own() that uses the linear kernel
#' and a unique variance for each environment to compute the full covariance
#' matrix and its partial derivatives.
#' @export
#'
mvlk <- function(C) {
  return(function(order, kappa) {
    # The correlation matrix of the traits:
    q <- nrow(C)
    p <- order / q
    Ct <- matrix(1, p, p)
    Ct[upper.tri(Ct)] <- Ct[lower.tri(Ct)] <- kappa[(order + 1):(order + ((p^2 - p) / 2))]

    # The full covariance matrix:
    S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
    V <- S * kronecker(Ct, C)

    # Derivatives wrt kappa[1:(p * q)] (variances):
    varderivs <- vector("list", order)
    for (dk in 1:order) {
      I <- matrix(0, order, order)
      I[dk,] <- I[, dk] <- 1
      tmp <- sqrt(kappa[1:order])
      tmp[dk] <- 1 / tmp[dk]
      tmp <- outer(tmp, tmp)
      tmp[dk, dk] <- 1
      deriv <- 0.5 * I * tmp * kronecker(Ct, C)
      deriv[dk, dk] <- 1
      varderivs[[dk]] <- deriv
    }

    # Derivatives wrt kappa[(p + 1):(p + ((p^2 - p) / 2))] (correlations):
    corderivs <- vector("list", (p^2 - p) / 2)
    for (dk in 1:((p^2 - p) / 2)) {
      # Indicator matrix of where kappa[dk] is present:
      I <- matrix(0, p, p)
      I[upper.tri(I)][dk] <- I[lower.tri(I)][dk] <- 1
      corderivs[[dk]] <- S * kronecker(I, C)
    }

    return(c(list(V), varderivs, corderivs))
  })
}
