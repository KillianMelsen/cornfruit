#' mvgk
#'
#' The \code{mvgk} function computes the covariance matrix and partial derivative matrices
#' for the multiple-variance Gaussian kernel model (MV-GK).
#' This function can be used for both GxE and GxExM/GxExT models.
#' The MV-GK model assumes that every environment by management combination has a
#' unique genetic variance. For GxE data, this means that the number of genetic
#' variances that are estimated is equal to the number of environments. Instead of
#' using a linear kernel, the MV-GK model uses a non-linear, Gaussian kernel for
#' which a bandwidth parameter is estimated.
#'
#' @param D A distance matrix in the form of a symmetric numeric matrix. Row-
#' and column-order must match the order of the levels of the environment
#' factor in the dataframe that will be passed to ASReml-R.
#'
#' @return A function compatible with asreml::own() that uses the distance matrix
#' to compute a non-linear kernel given a bandwidth \code{h}. This non-linear
#' kernel is then combined with unique variances for each environment to compute
#' the full covariance matrix and its partial derivatives.
#'
#' @examples
#' \dontrun{
#' # For GxE data:
#' vf <- cornfruit::mvgk(D = distanceMatrix)
#' q <- length(levels(df$Environment))
#' # Variances + bandwidth:
#' init <- c(rep(0.1, q), 0.1)
#' type <- c(rep("V", q), "V")
#' cons <- c(rep("P", q), "P")
#' fit <- asreml(fixed = Y ~ 1 + Environment,
#'               random = ~ own(Environment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExM data:
#' vf <- cornfruit::mvgk(D = distanceMatrix)
#' p <- length(levels(df$Management))
#' q <- length(levels(df$Environment))
#' # Variances + correlations + bandwidth:
#' init <- c(rep(0.1, p * q), rep(0.1, (p^2 - p) / 2), 0.1)
#' type <- c(rep("V", p * q), rep("R", (p^2 - p) / 2), "V")
#' cons <- c(rep("P", p * q), rep("U", (p^2 - p) / 2), "P")
#' fit <- asreml(fixed = Y ~ 1 + ManagementEnvironment,
#'               random = ~ own(ManagementEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExT data:
#' vf <- cornfruit::mvgk(D = distanceMatrix)
#' p <- length(levels(df$Trait))
#' q <- length(levels(df$Environment))
#' # Variances + correlations + bandwidth:
#' init <- c(rep(0.1, p * q), rep(0.1, (p^2 - p) / 2), 0.1)
#' type <- c(rep("V", p * q), rep("R", (p^2 - p) / 2), "V")
#' cons <- c(rep("P", p * q), rep("U", (p^2 - p) / 2), "P")
#' fit <- asreml(fixed = Y ~ 1 + TraitEnvironment,
#'               random = ~ own(TraitEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df,
#'               asmv = Trait)
#' }
#'
#' @export
#'
mvgk <- function(D) {
  return(function(order, kappa) {
    if (order == nrow(D)) {
      # if (all(kappa == get("init"))) {
      #   message("Fitting GxE MV-GK model...")
      # }
      # The full covariance matrix:
      S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
      C <- exp(-kappa[order + 1] * D)
      V <- S * C
      # Derivatives wrt kappa[1:order] (variances):
      varderivs <- vector("list", order)
      for (i in 1:order) {
        A <- matrix(0, order, order)
        A[, i] <- 0.5
        varderivs[[i]] <- V * (A + t(A)) / kappa[i]
      }
      # Covariance matrix, deriv wrt the variances and deriv wrt the bandwidth:
      return(c(list(V), varderivs, list(-D * V)))
    } else {
      # if (all(kappa == get("init"))) {
      #   message("Fitting GxExM MV-GK model...")
      # }
      # The correlation matrix of the traits:
      q <- nrow(D) # Number of environments
      p <- order / q # Number of managements
      Ct <- matrix(1, p, p)
      Ct[upper.tri(Ct)] <- Ct[lower.tri(Ct)] <- kappa[(order + 1):(order + ((p^2 - p) / 2))]

      # The full covariance matrix:
      S <- outer(sqrt(kappa[1:order]), sqrt(kappa[1:order]))
      C <- exp(-kappa[p * q + ((p^2 - p) / 2) + 1] * D)
      V <- S * kronecker(Ct, C)

      # Derivatives wrt kappa[1:(p * q)] (variances):
      varderivs <- vector("list", order)
      for (i in 1:order) {
        A <- matrix(0, order, order)
        A[, i] <- 0.5
        varderivs[[i]] <- V * (A + t(A)) / kappa[i]
      }

      # Derivatives wrt kappa[(p + 1):(p + ((p^2 - p) / 2))] (correlations):
      corderivs <- vector("list", (p^2 - p) / 2)
      for (i in 1:((p^2 - p) / 2)) {
        B <- matrix(0, p, p)
        B[upper.tri(B)][i] <- B[lower.tri(B)][i] <- 1
        corderivs[[i]] <- S * kronecker(B, C)
      }

      # Derivative wrt the bandwidth:
      dbw <- S * kronecker(Ct, -D * C)

      return(c(list(V), varderivs, corderivs, list(dbw)))
    }
  })
}
