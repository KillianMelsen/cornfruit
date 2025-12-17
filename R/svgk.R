#' svgk
#'
#' The \code{svgk} function computes the covariance matrix and partial derivative matrices
#' for the single-variance Gaussian kernel model (SV-GK).
#' This function can be used for both GxE and GxExM/GxExT models.
#' The SV-GK model assumes that all environments within a level of management or
#' trait have the same genetic variance. For GxE data, this means that only a
#' single genetic variance is estimated. Instead of using a linear kernel, the
#' SV-GK model uses a non-linear, Gaussian kernel for which a bandwidth parameter
#' is estimated.
#'
#' @param D A distance matrix in the form of a symmetric numeric matrix. Row-
#' and column-order must match the order of the levels of the environment
#' factor in the dataframe that will be passed to ASReml-R.
#'
#' @return A function compatible with asreml::own() that uses the distance matrix
#' to compute a non-linear, Gaussian kernel given a bandwidth \code{h}. This non-linear
#' kernel is then combined with a single variance for each environment to compute
#' the full covariance matrix and its partial derivatives.
#'
#' @examples
#' \dontrun{
#' # For GxE data:
#' vf <- cornfruit::svgk(D = distanceMatrix)
#' # Variance + bandwidth:
#' init <- c(0.1, 0.1)
#' type <- c("V", "V")
#' cons <- c("P", "P")
#' fit <- asreml(fixed = Y ~ 1 + Environment,
#'               random = ~ own(Environment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExM data:
#' vf <- cornfruit::svgk(D = distanceMatrix)
#' p <- length(levels(df$Management))
#' # Variances + correlations + bandwidth:
#' init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2), 0.1)
#' type <- c(rep("V", p), rep("R", (p^2 - p) / 2), "V")
#' cons <- c(rep("P", p), rep("U", (p^2 - p) / 2), "P")
#' fit <- asreml(fixed = Y ~ 1 + ManagementEnvironment,
#'               random = ~ own(ManagementEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExT data:
#' vf <- cornfruit::svgk(D = distanceMatrix)
#' p <- length(levels(df$Trait))
#' # Variances + correlations + bandwidth:
#' init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2), 0.1)
#' type <- c(rep("V", p), rep("R", (p^2 - p) / 2), "V")
#' cons <- c(rep("P", p), rep("U", (p^2 - p) / 2), "P")
#' fit <- asreml(fixed = Y ~ 1 + TraitEnvironment,
#'               random = ~ own(TraitEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df,
#'               asmv = Trait)
#' }
#'
#' @export
#'
svgk <- function(D) {
  return(function(order, kappa) {
    if (order == nrow(D)) {
      # if (all(kappa == get("init"))) {
      #   message("Fitting GxE SV-GK model...")
      # }
      # The full covariance matrix:
      C <- exp(-kappa[2] * D)
      V <- kappa[1] * C
      # Covariance matrix, deriv wrt the variance and deriv wrt the bandwidth:
      return(list(V, C, -D * V))
    } else {
      # if (all(kappa == get("init"))) {
      #   message("Fitting GxExM SV-GK model...")
      # }
      # The correlation matrix of the traits:
      q <- nrow(D) # Number of environments
      p <- order / q # Number of managements
      Ct <- matrix(1, p, p)
      Ct[upper.tri(Ct)] <- Ct[lower.tri(Ct)] <- kappa[(p + 1):(p + ((p^2 - p) / 2))]

      # The full covariance matrix:
      S <- outer(sqrt(kappa[1:p]), sqrt(kappa[1:p]))
      C <- exp(-kappa[p + ((p^2 - p) / 2) + 1] * D)
      V <- kronecker(S * Ct, C)

      # Derivatives wrt kappa[1:p] (variances):
      varderivs <- vector("list", p)
      for (i in 1:p) {
        A <- matrix(0, p, p)
        A[, i] <- 0.5
        varderivs[[i]] <- kronecker(S * Ct * (A + t(A)) / kappa[i], C)
      }

      # Derivatives wrt kappa[(p + 1):(p + ((p^2 - p) / 2))] (correlations):
      corderivs <- vector("list", (p^2 - p) / 2)
      for (i in 1:((p^2 - p) / 2)) {
        # Indicator matrix of where kappa[dk] is present:
        B <- matrix(0, p, p)
        B[upper.tri(B)][i] <- B[lower.tri(B)][i] <- 1
        corderivs[[i]] <- kronecker(S * B, C)
      }

      # Derivative wrt the bandwidth:
      bwderiv <- kronecker(S * Ct, -D * C)

      return(c(list(V), varderivs, corderivs, list(bwderiv)))
    }
  })
}
