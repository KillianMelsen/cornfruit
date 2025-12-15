#' DROPS GxExM data
#'
#' A list containing:
#' 1. A subset of the European Union's  DROPS (DROught-tolerant yielding PlantS) project dataset.
#' 2. A linear kernel of the included environments.
#' 3. A distance matrix of the included environments.
#' 4. A kinship matrix of the included genotypes.
#'
#' @format ## `exampleData$df`
#' A data frame with 300 rows and 5 columns:
#' \describe{
#'   \item{Man}{A factor for the management names (*R*ain-fed and *W*atered)}
#'   \item{Env}{A factor for the environment names (location-year)}
#'   \item{Gen}{A factor for the genotype names}
#'   \item{ManEnv}{A factor combining the levels of management and environment}
#'   \item{GY}{Grain Yield in ton per hectare}
#' }
#'
#' @format ## `exampleData$C`
#' A 3x3 linear kernel of the three included environments.
#'
#' @format ## `exampleData$D`
#' A 3x3 distance matrix of the three included environments to be used for the Gaussian kernels.
#'
#' @format ## `exampleData$K`
#' A 50x50 kinship matrix based on SNP markers for the 50 included genotypes.
#'
#' @source <https://doi.org/10.15454/IASSTN>
"exampleData"
