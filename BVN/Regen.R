set.seed(1234)

source("Initialization.R")
source("Functions.R")
load("Es.Rdata")

library(mcmcse)


norm_reg1 = matrix(0, nrow = reps, ncol = 4)
norm_reg2 = matrix(0, nrow = reps, ncol = 4)
norm_bmopt = matrix(0, nrow = reps, ncol = 4)
norm_bmth = matrix(0, nrow = reps, ncol = 4)

ess_reg1 = rep(0, reps)
ess_reg2 = rep(0, reps)
ess_bmopt = rep(0, reps)
ess_bmth = rep(0, reps)

for(k in 1:reps){
	dat = Gen_data(nsize)
	for(j in 1:length(samp_size)){
		work_d = dat[1:samp_size[j],]

		dummy1 = regen1_variance(work_d)
		norm_reg1[k,j] = dummy1$Fnorm

		dummy2 = regen2_variance(work_d)
		norm_reg2[k,j] = dummy2$Fnorm

		dummy3 = mcse.multi(work_d[,-3], method = "bm", r = 1)$cov
		norm_bmopt[k,j] = norm(Tr - dummy3, type = "F")

		dummy4 = mcse.multi(work_d[,-3], method = "bm", size = floor((1e4)^(1/2 + .00001)),
		 r = 1)$cov
		norm_bmth[k,j] = norm(Tr - dummy4, type = "F")

		if(samp_size[j] == max(samp_size)){
			ess_reg1[k] = dummy1$ESS
			ess_reg2[k] = dummy2$ESS
			ess_bmopt[k] = multiESS(work_d[,-3], covmat = dummy3)
			ess_bmth[k] = multiESS(work_d[,-3], covmat = dummy4)
		}
		print(samp_size[j])
	}
	print(k)
}

save(norm_reg1, ess_reg1, norm_reg2, ess_reg2, norm_bmopt, ess_bmopt, norm_bmth, ess_bmth, file = "NE.Rdata")








