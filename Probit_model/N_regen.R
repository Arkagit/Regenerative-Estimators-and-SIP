set.seed(4321)
source("Data_gen.R")
source("Small_set.R")


# Generating samples from AC algorithm
nsim <- 1e6
ac.out.1e6 <- probit_gibbs(dat = dat, nsim = nsim)

#dim(ac.out.1e6$beta)
#Initialization for eta prob
eta1_det = numeric(length = nsim)
eta1_randm = numeric(length = nsim)
# Cuboidal regions
for(j in 1:(nsim-1)){
  if(j %% (nsim/10) == 0) print(paste(j/nsim * 100, "% done"))
  betaj1 <- ac.out.1e6$beta[j+1,]
  zj1 <- ac.out.1e6$z[j+1,]
  betaj0 <- ac.out.1e6$beta[j,]
  zj0 <- ac.out.1e6$z[j,]
  eta1_det[j] <- etaj_cube_det(betaj0, zj0, betaj1, zj1)
  eta1_randm[j] <- etaj_cube_randm(betaj0, zj0, betaj1, zj1)
}

#Number of total regenerations for Deterministic Scan
sum_regen_det = rbinom(nsim, 1, eta1_det)
sum(sum_regen_det)

#Number of total regenerations for Random Scan
sum_regen_randm = rbinom(nsim, 1, eta1_randm)
sum(sum_regen_randm)