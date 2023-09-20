set.seed(4321)
library(mcmcse)
library(MASS)
library(truncnorm)

source("Data_gen.R")
source("probit_chain.R")


init <- 1e6
counter <- numeric(0)

for(i in 1:100){
	out.1e6 <- probit_gibbs(dat = dat, nsim = init)
	counter[i] = batchSize(out.1e6$beta, method = "bm")/(init^(1/3))
}

ratio = floor(mean(counter)); ratio