set.seed(1234)
source("Data_gen.R")
source("probit_chain.R")
load("Small.Rdata")


# Generating samples from AC algorithm
nsim <- 1e6
ac.out.1e6 <- probit_gibbs(dat = dat, nsim = nsim)

#dim(ac.out.1e6$beta)
#Initialization for eta prob
eta1_det = numeric(length = nsim)
eta2_randm = numeric(length = nsim)
# Cuboidal regions
for(j in 1:(nsim-1)){
  if(j %% (nsim/10) == 0) print(paste(j/nsim * 100, "% done"))
  betaj1 <- ac.out.1e6$beta[j+1,]
  zj1 <- ac.out.1e6$z[j+1,]
  betaj0 <- ac.out.1e6$beta[j,]
  zj0 <- ac.out.1e6$z[j,]
  eta1_det[j] <- etaj_cube_det(betaj0, zj0, betaj1, zj1)
}

#Number of total regenerations for Deterministic Scan
sum_regen_det = rbinom(nsim, 1, eta1_det)
sum(sum_regen_det)
hist(eta1_det)


for(j in 1:(nsim-2)){
  if(j %% (nsim/10) == 0) print(paste(j/nsim * 100, "% done"))
  betaj2 <- ac.out.1e6$beta[j+2,]
  zj2 <- ac.out.1e6$z[j+2,]
  betaj0 <- ac.out.1e6$beta[j,]
  zj0 <- ac.out.1e6$z[j,]
  eta2_randm[j] <- etaj_cube_randm2(betaj0, zj0, betaj2, zj2)
}


#Number of total regenerations for Random Scan
sum_regen_randm = rbinom(nsim, 1, eta2_randm)
sum(sum_regen_randm)