set.seed(1234)

library(mcmcse)
library(foreach)
library(doParallel)

load("Small.Rdata")
source("Data_gen.R")
source("probit_chain.R")

reps = 10
Ratio = 15
samp_size = 1e6
# Parallelizing norm calculation

cv = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

#register it to be used by %dopar%
doParallel::registerDoParallel(cores = n.cores)
k=1

covar = foreach(k = 1:reps, .packages = c("mcmcse","truncnorm")) %dopar% {
#covar  = list()


#for(k in 1:reps){

	ac_data = probit_gibbs(dat = dat, nsim = samp_size)
	eta1_det = rep(0, samp_size)

	for(j in 1:(samp_size-1)){
	  betaj1 <- ac_data$beta[j+1,]
	  zj1 <- ac_data$z[j+1,]
	  betaj0 <- ac_data$beta[j,]
	  zj0 <- ac_data$z[j,]
	  eta1_det[j] <- etaj_cube_det(betaj0, zj0, betaj1, zj1)
	}
	#Number of total regenerations for Deterministic Scan
	regen_det = rbinom(samp_size, 1, eta1_det)

	# Summary on regenerations
	regen_steps = which(regen_det == 1) + 1
	regen_points = c(1, regen_steps[- length(regen_steps)])
	regen_length = regen_steps - regen_points

	Z = matrix(0, nrow = length(regen_steps), ncol = dim(x)[2])

	for (j in 1:length(regen_steps)) {
	    foo <- ac_data$beta[regen_points[j]:(regen_points[j] + regen_length[j]), ]
	    foo <- as.matrix(foo, ncol = 3)
	    Z[j,] = colSums(foo)
	}

	Z_bar = colMeans(Z)

	beta_reg_mean = colMeans(ac_data$beta[1:max(regen_steps),])

	# Mean length of regenerations
	mu = mean(regen_length)

	dummy = list()

	dummy[[1]] = mcse.multi(ac_data$beta, method = "bm", r = 1, size = "sqroot")$cov

	dummy[[2]] = mcse.multi(ac_data$beta, method = "bm", r = 1)$cov

	dummy[[3]] = mcse.multi(ac_data$beta, method = "bm", r = 1, size = 5e4)$cov

	dummy[[4]] = regen_var(Z, mu, regen_steps)

	dummy[[5]] = regen_var2(Z, mu, regen_steps)

	dummy[[6]] = regen_var3(Z, mu, regen_steps)

	cv[[k]] = dummy

	cv[[k]]

	#print(paste(k))
}

covar


norm_reg1 = rep(0, length(reps))
norm_reg2 = rep(0, length(reps))
norm_reg3 = rep(0, length(reps))
norm_bmopt = rep(0, length(reps))
norm_bmth = rep(0, length(reps))
norm_bm1e4 = rep(0, length(reps))


for(i in 1:reps){
	norm_bmth[i] = norm(cv[[i]][[1]], type = "F")
	norm_bmopt[i] = norm(cv[[i]][[2]], type = "F")
	norm_bm1e4[i] = norm(cv[[i]][[3]], type = "F")
	norm_reg1[i] = norm(cv[[i]][[4]], type = "F")
	norm_reg2[i] = norm(cv[[i]][[5]], type = "F")
	norm_reg3[i] = norm(cv[[i]][[6]], type = "F")
}

names = c("BM(0.5)", "BM(1/3)", "BM(1e4)", "REG1", "REG2", "REG3")


pdf(file="Boxplot1.pdf")

boxplot(norm_bmth, norm_bmopt, norm_bm1e4, norm_reg1, norm_reg2, norm_reg3,
	names = names)

dev.off()



