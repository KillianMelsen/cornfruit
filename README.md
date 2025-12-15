### Installation
```
devtools::install_github("KillianMelsen/cornfruit")
library(cornfruit)
library(asreml)

df <- exampleData$df
K <- exampleData$K
D <- exampleData$D
p <- length(levels(df$Man))
q <- length(levels(df$Env))

vf <- mvgk(D = D)
init <- (rep(0.1, p * q), rep(0.1, (p^2 - p) / 2), 0.1)
type <- (rep("V", p * q), rep("R", (p^2 - p) / 2), "V")
cons <- (rep("P", p * q), rep("U", (p^2 - p) / 2), "P")

fit <- asreml(fixed = GY ~ -1 + ManEnv,
              random = ~ own(ManEnv, "vf", init, type, cons):vm(Gen, K),
              residual = ~ units,
              data = df)
```
