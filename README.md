### Installation
```R
devtools::install_github("KillianMelsen/cornfruit")
```
### Usage
```R
# Packages:
library(cornfruit)
library(asreml)

# Loading data and determining number of managemens and environments:
df <- exampleData$df # Data frame
K <- exampleData$K # kinship
disMatrix <- exampleData$D # Distance Matrix for Gaussian kernel
linKernel <- exampleData$C # Linear kernel
p <- length(levels(df$Man)) # Number of managements
q <- length(levels(df$Env)) # Number of environments

# GxExM models =================================================================
## SV-LK model =================================================================
# Defining the variance function:
vf <- svlk(C = linKernel)

# Setting initial values, parameter types, and constraints.
# SV-LK, so 1 variance for each of the p managements, and (p^2 - p) / 2
# correlations between the p managements:
init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2))
type <- c(rep("V", p), rep("R", (p^2 - p) / 2))
cons <- c(rep("P", p), rep("U", (p^2 - p) / 2))

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## SV-GK model =================================================================
# Defining the variance function:
vf <- svgk(D = disMatrix)

# Setting initial values, parameter types, and constraints.
# SV-GK, so 1 variance for each of the p managements, (p^2 - p) / 2 correlations
# between the p managements, and a single bandwidth:
init <- c(rep(0.1, p), rep(0.1, (p^2 - p) / 2), 0.1)
type <- c(rep("V", p), rep("R", (p^2 - p) / 2), "V")
cons <- c(rep("P", p), rep("U", (p^2 - p) / 2), "P")

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## MV-LK model =================================================================
# Defining the variance function:
vf <- mvlk(C = linKernel)

# Setting initial values, parameter types, and constraints.
# MV-LK, so 1 variance for each of the p * q management by environment
# combinations, and (p^2 - p) / 2 correlations between the p managements:
init <- c(rep(0.1, p * q), rep(0.1, (p^2 - p) / 2))
type <- c(rep("V", p * q), rep("R", (p^2 - p) / 2))
cons <- c(rep("P", p * q), rep("U", (p^2 - p) / 2))

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## MV-GK model =================================================================
# Defining the variance function:
vf <- mvgk(D = disMatrix)

# Setting initial values, parameter types, and constraints.
# MV-GK, so 1 variance for each of the p * q management by environment
# combinations, (p^2 - p) / 2 correlations between the p managements, and a
# single bandwidth:
init <- c(rep(0.1, p * q), rep(0.1, (p^2 - p) / 2), 0.1)
type <- c(rep("V", p * q), rep("R", (p^2 - p) / 2), "V")
cons <- c(rep("P", p * q), rep("U", (p^2 - p) / 2), "P")

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

# GxE models ===================================================================
# Dropping one management:
df <- droplevels(df[df$Man == "R",])

## SV-LK model =================================================================
# Defining the variance function:
vf <- svlk(C = linKernel)

# Setting initial values, parameter types, and constraints.
# SV-LK, so 1 variance for all environments:
init <- 0.1
type <- "V"
cons <- "P"

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ own(Env, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## SV-GK model =================================================================
# Defining the variance function:
vf <- svgk(D = disMatrix)

# Setting initial values, parameter types, and constraints.
# SV-GK, so 1 variance for all environments, and a single bandwidth:
init <- c(0.1, 0.1)
type <- c("V", "V")
cons <- c("P", "P")

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ own(Env, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## MV-LK model =================================================================
# Defining the variance function:
vf <- mvlk(C = linKernel)

# Setting initial values, parameter types, and constraints.
# MV-LK, so 1 variance for each of the q environments:
init <- rep(0.1, q)
type <- rep("V", q)
cons <- rep("P", q)

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ own(Env, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)

## MV-GK model =================================================================
# Defining the variance function:
vf <- mvgk(D = disMatrix)

# Setting initial values, parameter types, and constraints.
# MV-GK, so 1 variance for each of the q environments, and a single bandwidth:
init <- c(rep(0.1, q), 0.1)
type <- c(rep("V", q), "V")
cons <- c(rep("P", q), "P")

# Fitting the model:
fit <- asreml(fixed = GY ~ -1 + Env,
              random = ~ own(Env, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)
```
