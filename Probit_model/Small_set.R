set.seed(4321)

library(MASS)
library(truncnorm)
source("Data_gen.R")
source("probit_chain.R")

###############################################################################
###############################################################################
## Calculation of D* 

init <- 1e6
out.1e6 <- probit_gibbs(dat = dat, nsim = init)

tuner <- 0.09
z.star <- colMeans(out.1e6$z[1:init, ])
beta.bar <- colMeans(out.1e6$beta[1:init, ])
s <- apply(out.1e6$beta[1:init, ], 2, sd)
c <- beta.bar - tuner*s
d <- beta.bar + tuner*s
V <- var(out.1e6$beta)
Vinv <- solve(V)




