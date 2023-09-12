library(mcmcse)
library(foreach)
library(doParallel)


source("Initialization.R")
source("Functions.R")

bs = numeric()

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

bs = foreach(i = 1:100, .packages = c("mcmcse"), .combine = "c")%dopar%{
	dat = Gen_data(1e6)
	a = batchSize(dat[,-3], method = "bm")
}

bs

const = ceiling(mean(bs))/(1e6)^(1/3); const

