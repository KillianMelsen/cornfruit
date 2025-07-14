#' mvgk
#'
#' @param D A distance matrix in the form of a symmetric numeric matrix. Row-
#' and column-order must match the order of the levels of the environment
#' factor in the dataframe that will be passed to ASReml-R.
#'
#' @return A function compatible with asreml::own() that uses the distance matrix
#' to compute a non-linear kernel given a bandwidth \code{h}. This non-linear
#' kernel is then combined with unique variances for each environment to compute
#' the full covariance matrix and its partial derivatives.
#' @export
#'
mvgk <- function(D) {
  return(function(order, kappa) {
    # The correlation matrix of the traits:
    q <- nrow(D)
    p <- order / q
    Ct <- matrix(1, p, p)
    Ct[upper.tri(Ct)] <- Ct[lower.tri(Ct)] <- kappa[(order + 1):(order + ((p^2 - p) / 2))]

    # The full covariance matrix:
    S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
    V <- S * kronecker(Ct, exp(-kappa[p * q + ((p^2 - p) / 2) + 1] * D))

    # Derivatives wrt kappa[1:(p * q)] (variances):
    varderivs <- vector("list", order)
    for (dk in 1:order) {
      I <- matrix(0, order, order)
      I[dk,] <- I[, dk] <- 1
      tmp <- sqrt(kappa[1:order])
      tmp[dk] <- 1 / tmp[dk]
      tmp <- outer(tmp, tmp)
      tmp[dk, dk] <- 1
      deriv <- 0.5 * I * tmp * kronecker(Ct, exp(-kappa[p * q + ((p^2 - p) / 2) + 1] * D))
      deriv[dk, dk] <- 1
      varderivs[[dk]] <- deriv
    }

    # Derivatives wrt kappa[(p + 1):(p + ((p^2 - p) / 2))] (correlations):
    corderivs <- vector("list", (p^2 - p) / 2)
    for (dk in 1:((p^2 - p) / 2)) {
      # Indicator matrix of where kappa[dk] is present:
      I <- matrix(0, p, p)
      I[upper.tri(I)][dk] <- I[lower.tri(I)][dk] <- 1
      corderivs[[dk]] <- S * kronecker(I, exp(-kappa[p * q + ((p^2 - p) / 2) + 1] * D))
    }

    # Derivative wrt the bandwidth:
    dbw <- S * kronecker(Ct, -D * exp(-kappa[p * q + ((p^2 - p) / 2) + 1] * D))

    return(c(list(V), varderivs, corderivs, list(dbw)))
  })
}
