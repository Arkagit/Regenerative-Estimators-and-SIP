set.seed(1111)

source("Initialization.R")
source("Functions.R")
load("Es.Rdata")

library(mcmcse)
library(foreach)
library(doParallel)

reps = 1000
const = 0.25
# Parallelizing norm calculation

cv = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

#my.cluster <- parallel::makeCluster(
#  n.cores, 
#  type = "PSOCK"
#  )

#check cluster definition (optional)
#print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cores = n.cores)

#check if it is registered (optional)
#foreach::getDoParRegistered()

# rep -> samp_size -> estimator
#covar = list()

covar = foreach(k = 1:reps, .packages = c("mcmcse")) %dopar% {
#for(k in 1:reps){
	dat = Gen_data(max(samp_size))
	regens2 = rbinom(max(samp_size), 1, Es * dat[,3])
	regens1 = rbinom(max(samp_size), 1, dat[,3])


	for(j in 1:length(samp_size)){
		work_d = dat[1:samp_size[j],]
		work_regen1 = regens1[1:samp_size[j]]
		work_regen2 = regens2[1:samp_size[j]]

		dummy = list()

		dummy[[1]] = regen1_variance(work_d, work_regen1)$cov

		dummy[[2]] = regen2_variance(work_d, work_regen2)$cov

		dummy[[3]] = mcse.multi(work_d[,-3], method = "bm", r = 1)$cov

		dummy[[4]] = mcse.multi(work_d[,-3], method = "bm", r = 1, size = floor(const*(samp_size[j])^(1/2 + .00001)))$cov

		cv[[j]] = dummy
	}
	cv
	#print(paste(k))
}

covar

norm_reg1 = matrix(0, nrow = reps, ncol = length(samp_size))
norm_reg2 = matrix(0, nrow = reps, ncol = length(samp_size))
norm_bmopt = matrix(0, nrow = reps, ncol = length(samp_size))
norm_bmth = matrix(0, nrow = reps, ncol = length(samp_size))


for(i in 1:reps){
	for(j in 1:length(samp_size)){
		norm_reg1[i,j] = norm(Tr - covar[[i]][[j]][[1]], type = "F")
		norm_reg2[i,j] = norm(Tr - covar[[i]][[j]][[2]], type = "F")
		norm_bmopt[i,j] = norm(Tr - covar[[i]][[j]][[3]], type = "F")
		norm_bmth[i,j] = norm(Tr - covar[[i]][[j]][[4]], type = "F")
	}
	
}

# Parallelizing ESS calculation

cv = list()

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

#my.cluster <- parallel::makeCluster(
#  n.cores, 
#  type = "PSOCK"
#  )

#check cluster definition (optional)
#print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cores = n.cores)

#check if it is registered (optional)
#foreach::getDoParRegistered()
# rep -> samp_size -> estimator

#ESS = list()

ESS = foreach(k = 1:reps, .packages = c("mcmcse")) %dopar% {
#for(k in 1:reps){
	dat = Gen_data(max(samp_size))
	regens2 = rbinom(max(samp_size), 1, Es * dat[,3])
	regens1 = rbinom(max(samp_size), 1, dat[,3])
	for(j in 1:length(samp_size)){
		work_d = dat[1:samp_size[j],]; dim(work_d)
		work_regen1 = regens1[1:samp_size[j]]
		work_regen2 = regens2[1:samp_size[j]]

		dummy = list()

		dummy[[1]] = multiESS(work_d[,-3], covmat = regen1_variance(work_d, work_regen1)$cov)/samp_size[j]

		dummy[[2]] = multiESS(work_d[,-3], covmat = regen2_variance(work_d, work_regen2)$cov)/samp_size[j]
		
		dummy[[3]] = multiESS(work_d[,-3], covmat = mcse.multi(work_d[,-3], method = "bm", r = 1)$cov)/samp_size[j]
		
		dummy[[4]] = multiESS(work_d[,-3], covmat = mcse.multi(work_d[,-3], method = "bm", r = 1, size = floor(const*(samp_size[j])^(1/2 + .00001)))$cov)/samp_size[j]

		cv[[j]] = dummy
	}
	cv
	#print(paste(k))
}


ESS


ess_reg1 = matrix(0, nrow = reps, ncol = length(samp_size))
ess_reg2 = matrix(0, nrow = reps, ncol = length(samp_size))
ess_bmopt = matrix(0, nrow = reps, ncol = length(samp_size))
ess_bmth = matrix(0, nrow = reps, ncol = length(samp_size))


for(i in 1:reps){
	for(j in 1:length(samp_size)){
		ess_reg1[i,j] = ESS[[i]][[j]][[1]]
		ess_reg2[i,j] = ESS[[i]][[j]][[2]]
		ess_bmopt[i,j] = ESS[[i]][[j]][[3]]
		ess_bmth[i,j] = ESS[[i]][[j]][[4]]
	}
	
}

save(reps, norm_reg1, ess_reg1, norm_reg2, ess_reg2, norm_bmopt, ess_bmopt, norm_bmth, ess_bmth, file = "NE.Rdata")








