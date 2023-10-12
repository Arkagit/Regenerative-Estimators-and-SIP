set.seed(1234)

library(mcmcse)
library(truncnorm)
library(foreach)
library(doParallel)

load("Small.Rdata")
source("Data_gen.R")
source("probit_chain.R")

reps = 1000
#Ratio = 15
samp_size = c(1e4, 5e4, 1e5, 3e5, 5e5, 8e5, 1e6)
# Parallelizing norm calculation

cv = list()

#parallel::detectCores()
#n.cores <- parallel::detectCores() - 1

#register it to be used by %dopar%
#doParallel::registerDoParallel(cores = n.cores)

EST  = list()


for(k in 1:reps){
#covar = foreach(k = 1:reps, .packages = c("mcmcse","truncnorm")) %dopar% {

	ac_data = probit_gibbs(dat = dat, nsim = max(samp_size))

	foo = list()

	for(i in 1:length(samp_size)){

		foo_data = ac_data$beta[1:samp_size[i],]

		dummy = list()

		dummy[[1]] = mcse.multi(foo_data, method = "bm", r = 1, size = (samp_size[i])^(1/2 + 0.0001))$cov

		dummy[[2]] = mcse.multi(foo_data, method = "bm", r = 1)$cov

		dummy[[3]] = multiESS(foo_data, dummy[[1]])

		dummy[[4]] = multiESS(foo_data, dummy[[2]])

		foo[[i]] = dummy

		print(paste(i))
	}

	#foo_ls
	EST[[k]] = foo
	
	print(paste(k))
}

EST

#covar = covar[-which(sapply(covar, is.null))]

###############################

#ESS = list()

#

#for(k in 1:reps){
#ESS = foreach(k = 1:reps, .packages = c("mcmcse","truncnorm")) %dopar% {

#	ac_data = probit_gibbs(dat = dat, nsim = max(samp_size))
#
#	foo_es = list()

#	for(i in 1:length(samp_size)){

#		foo_data = ac_data$beta[1:samp_size[i],]

#		dummy = list()
#
#		dummy[[1]] = multiESS(foo_data, mcse.multi(foo_data, method = "bm", r = 1, size = (samp_size[i])^(1/2 + 0.0001))$cov)
#
#		dummy[[2]] = multiESS(foo_data, mcse.multi(foo_data, method = "bm", r = 1)$cov)

#		foo_es[[i]] = dummy
#
#		print(paste(i))
#	}

	#foo_es
#	ESS[[k]] = foo_es
#	print(paste(k))
#}

#ESS

#ESS = ESS[-which(sapply(ESS, is.null))]

save(EST, samp_size, reps, file = "variance.Rdata")




#names = c("BM(1/2 + 0.0001)", "BM(1/3)")#, "BM(1e4)", "REG1", "REG2", "REG3")


#pdf(file="Boxplot1.pdf")

#boxplot(norm_bmth, norm_bmopt, norm_bm1e4, norm_reg1, norm_reg2, norm_reg3,
#	names = names)

#dev.off()



