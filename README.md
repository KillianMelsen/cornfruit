### Installation
```R
devtools::install_github("KillianMelsen/cornfruit")
```
### Usage
```R
library(cornfruit)
library(asreml)

df <- exampleData$df
K <- exampleData$K
D <- exampleData$D
p <- length(levels(df$Man))
q <- length(levels(df$Env))

vf <- mvgk(D = D)
init <- c(rep(0.1, p * q), rep(0.1, (p^2 - p) / 2), 0.1)
type <- c(rep("V", p * q), rep("R", (p^2 - p) / 2), "V")
cons <- c(rep("P", p * q), rep("U", (p^2 - p) / 2), "P")

fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)
```
