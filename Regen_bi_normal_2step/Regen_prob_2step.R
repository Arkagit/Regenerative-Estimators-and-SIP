set.seed(1234)

source("Mean_calc.R")
library(mcmcse)

# Initialization
nsim = 1e6
x = rep(0, nsim)
y = rep(0, nsim)
s = rep(0, nsim)


x[1] = 5
y[1] = rnorm(1, 2*x[1], 1)
s[1] = (pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))) * exp(-(4*y[1]^{2} + 100*y[1])/100)


for(i in 2:nsim){
	x[i] = rnorm(1, 2*y[i-1]/5, 1/sqrt(5))
	y[i] = rnorm(1, 2*x[i], 1)
	if(y[i] > 0){
		s[i] = (pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))) * exp(-(2*y[i]^{2} + 50*y[i])/5)
	}else{
		s[i] = (pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))) * exp(-(2*y[i]^{2} - 50*y[i])/5)
	}	
}

new_s = Es*s

# Minorization constants/Regeneration probabilities
hist(new_s)

# Regeneration counters
regens = rbinom(nsim, 1, new_s)

# Total number of detected regenerations
sum(regens)

T = append(c(0),which(regens == 1)); head(T)
Z = matrix(0, nrow = sum(regens), ncol = 2); head(Z)

for(i in 1:(length(T)-1)){
	for(j in (T[i] + 1):T[i+1]){
		Z[i,1] = Z[i,1] + x[j]
		Z[i,2] = Z[i,2] + y[j]
	}
}

# Individual regeneration times
regen_times = (which(regens == 1)) - append(c(0), which(regens == 1)[-length(which(regens == 1))])

# Comparing regeneration mean and sample mean
colSums(Z)/sum(regen_times); c(mean(x), mean(y))

# Construction of function for regeneration
variance = function(Z, regen_times){
	R = length(regen_times)
	Z_bar = colMeans(Z)
	A1 = matrix(0, nrow= 2, ncol = 2)
	A2 = matrix(0, nrow= 2, ncol = 2)
	A3 = matrix(0, nrow= 2, ncol = 2)
	for(i in 1:R){
		A1 = A1 + (Z[i,] - Z_bar)%*%t(Z[i,] - Z_bar)
	}
	for (i in 1:(R-1)) {
		A2 = A2 + (Z[i,] - Z_bar)%*%t(Z[i+1,] - Z_bar)
		A3 = A3 + (Z[i+1,] - Z_bar)%*%t(Z[i,] - Z_bar)
	}
	V = A1 + A2 + A3
	return(V/sum(regen_times))
}


# Dataset
dataset = cbind(x, y); head(dataset)

# Comparing regeneration variance with other methods
variance(Z, regen_times)
mcse.multi(dataset, method = "bm", r = 1)$cov
mcse.multi(dataset, method = "bm", r = 1, size = floor(nsim^(1/2 + .00001)))$cov
mcse.initseq(dataset)$cov
