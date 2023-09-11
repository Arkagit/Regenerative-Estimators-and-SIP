set.seed(1234)

source("Initialization.R")
source("Functions.R")
load("Es.Rdata")

library(mcmcse)
library(foreach)
library(doParallel)

norm_reg1 = rep(0, length(samp_size))
norm_reg2 = rep(0, length(samp_size))
norm_bmopt = rep(0, length(samp_size))
norm_bmth = rep(0, length(samp_size))

ess_reg1 = rep(0, reps)
ess_reg2 = rep(0, reps)
ess_bmopt = rep(0, reps)
ess_bmth = rep(0, reps)

reps = 10

parallel::detectCores()
n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
  )

#check cluster definition (optional)
print(my.cluster)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

norm = foreach(k = 1:reps) %do% {
	dat = Gen_data(max(samp_size))
	for(j in 1:length(samp_size)){
		work_d = dat[1:samp_size[j],]

		dummy1 = regen1_variance(work_d)
		norm_reg1[j] = dummy1$Fnorm

		dummy2 = regen2_variance(work_d)
		norm_reg2[j] = dummy2$Fnorm

		dummy3 = mcse.multi(work_d[,-3], method = "bm", r = 1)$cov
		norm_bmopt[j] = norm(Tr - dummy3, type = "F")

		dummy4 = mcse.multi(work_d[,-3], method = "bm", size = floor((samp_size[j])^(1/2 + .00001)),
		 r = 1)$cov
		norm_bmth[j] = norm(Tr - dummy4, type = "F")

		#if(samp_size[j] == max(samp_size)){
			#ess_reg1[k] = dummy1$ESS
			#ess_reg2[k] = dummy2$ESS
			#ess_bmopt[k] = multiESS(work_d[,-3], covmat = dummy3)
			#ess_bmth[k] = multiESS(work_d[,-3], covmat = dummy4)
		#}
		#print(samp_size[j])
	}
	rbind(norm_reg1, norm_reg2, norm_bmopt, norm_bmth)

	#print(k)
}

norm

norm_reg1 = matrix(0, nrow = reps, ncol = 4)
norm_reg2 = matrix(0, nrow = reps, ncol = 4)
norm_bmopt = matrix(0, nrow = reps, ncol = 4)
norm_bmth = matrix(0, nrow = reps, ncol = 4)

for(i in 1:reps){
	norm_reg1[i,] = norm[[i]][1,1:length(samp_size)]
	norm_reg2[i,] = norm[[i]][2,1:length(samp_size)]
	norm_bmopt[i,] = norm[[i]][3,1:length(samp_size)]
	norm_bmth[i,] = norm[[i]][4,1:length(samp_size)]
}

for(k in 1:reps){
	dat = Gen_data(nsize)

	ess_reg1[k] = regen1_variance(dat)$ESS
	ess_reg2[k] = regen2_variance(dat)$ESS

	dummy3 = mcse.multi(dat[,-3], method = "bm", r = 1)$cov
	ess_reg1[k] = multiESS(dat[,-3], covmat = dummy3)

	dummy4 = mcse.multi(dat[,-3], method = "bm", size = floor((nsize)^(1/2 + .00001)),
		 r = 1)$cov
	ess_reg1[k] = multiESS(dat[,-3], covmat = dummy4)
}

save(norm_reg1, ess_reg1, norm_reg2, ess_reg2, norm_bmopt, ess_bmopt, norm_bmth, ess_bmth, file = "NE.Rdata")








