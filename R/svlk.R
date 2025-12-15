#' svlk
#'
#' The \code{svlk} function computes the covariance matrix and partial derivative matrices
#' for the single-variance linear kernel model (SV-LK).
#' This function can be used for both GxE and GxExM/GxExT models.
#' The SV-LK model assumes that all environments within a level of management or
#' trait have the same genetic variance. For GxE data, this means that only a
#' single genetic variance is estimated.
#'
#' @param C A linear kernel in the form of a symmetric numeric matrix. Row-
#' and column-order must match the order of the levels of the environment
#' factor in the dataframe that will be passed to ASReml-R.
#'
#' @return A function compatible with asreml::own() that uses the linear kernel
#' and a single variance for all environments (within a level of management) to
#' compute the full covariance matrix and its partial derivatives.
#'
#' @examples
#' \dontrun{
#' # For GxE data:
#' vf <- cornfruit::svlk(C = linearKernel)
#' # Variance only:
#' init <- 0.1
#' type <- "V"
#' cons <- "P"
#' fit <- asreml(fixed = Y ~ 1 + Environment,
#'               random = ~ own(Environment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExM data:
#' vf <- cornfruit::svlk(C = linearKernel)
#' p <- length(levels(df$Management))
#' # Variances + correlations:
#' init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2))
#' type <- c(rep("V", p), rep("R", (p^2 - p) / 2))
#' cons <- c(rep("P", p), rep("U", (p^2 - p) / 2))
#' fit <- asreml(fixed = Y ~ 1 + ManagementEnvironment,
#'               random = ~ own(ManagementEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df)
#'
#' # For GxExT data:
#' vf <- cornfruit::svlk(C = linearKernel)
#' p <- length(levels(df$Trait))
#' # Variances + correlations:
#' init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2))
#' type <- c(rep("V", p), rep("R", (p^2 - p) / 2))
#' cons <- c(rep("P", p), rep("U", (p^2 - p) / 2))
#' fit <- asreml(fixed = Y ~ 1 + TraitEnvironment,
#'               random = ~ own(TraitEnvironment, "vf", init, type, cons):vm(Genotype, K),
#'               residual = ~ units,
#'               data = df,
#'               asmv = Trait)
#' }
#'
#' @export
#'
svlk <- function(C) {
  return(function(order, kappa) {
    if (order == nrow(C)) {
      # The full covariance matrix:
      V <- kappa[1] * C
      return(list(V, C))
    } else {
      # The correlation matrix of the traits/managements:
      q <- nrow(C) # Number of environments
      p <- order / q # Number of managements
      Ct <- matrix(1, p, p)
      Ct[upper.tri(Ct)] <- Ct[lower.tri(Ct)] <- kappa[(p + 1):(p + ((p^2 - p) / 2))]

      # The full covariance matrix:
      S <- outer(sqrt(kappa[1:p]), sqrt(kappa[1:p]))
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
        B <- matrix(0, p, p)
        B[upper.tri(B)][i] <- B[lower.tri(B)][i] <- 1
        corderivs[[i]] <- kronecker(S * B, C)
      }

      return(c(list(V), varderivs, corderivs))
    }
  })
}
