set.seed(1234)

library(mcmcse)

source("Initialization.R")
load("Es.Rdata")

variance = function(Z, regen_times){
	R = length(regen_times)
	tau_bar = mean(regen_times)
	Z_bar = colMeans(Z)
	A1 = matrix(0, nrow= 2, ncol = 2)
	A2 = matrix(0, nrow= 2, ncol = 2)
	A3 = matrix(0, nrow= 2, ncol = 2)
	for(i in 1:R){
		A1 = A1 + (Z[i,] - Z_bar*regen_times[i]/tau_bar)%*%t(Z[i,] - Z_bar*regen_times[i]/tau_bar)
	}
	for (i in 1:(R-1)) {
		A2 = A2 + (Z[i,] - Z_bar*regen_times[i]/tau_bar)%*%t(Z[i+1,] - Z_bar*regen_times[i+1]/tau_bar)
		A3 = A3 + (Z[i+1,] - Z_bar*regen_times[i+1]/tau_bar)%*%t(Z[i,] - Z_bar*regen_times[i]/tau_bar)
	}
	V = A1 + A2 + A3
	return(V/sum(regen_times))
}

regen2_variance = function(xydata, regens){
	#l = dim(xydata)[1]
	#regens = rbinom(l, 1, Es * xydata[,3])

	#save(dat, "bug_test1.Rdata")

	# Individual regeneration times
	regen_times = (which(regens == 1)) -
	         append(c(0), which(regens == 1)[-length(which(regens == 1))])

	mu = mean(regen_times)

	T = append(c(0),which(regens == 1))
	Z = matrix(0, nrow = sum(regens), ncol = 2)

	for(i in 1:(length(T)-1)){
		for(j in (T[i] + 1):T[i+1]){
			Z[i,1] = Z[i,1] + xydata[j,1]
			Z[i,2] = Z[i,2] + xydata[j,2]
		}
	}

	Z_bar = colMeans(Z)

	W = Z - (regen_times %*% t(Z_bar))/mu

	#ess = multiESS(xydata[,-3], covmat = variance(Z, regen_times))
	#regVAR2 = norm(Tr - variance(Z, regen_times), type = "F")
	cov = variance(W, regen_times)
	return(list("cov" = cov))
}

regen1_variance = function(xydata, regens){
	#l = dim(xydata)[1]
	#regens = rbinom(l, 1, xydata[,3])

	#save(dat, "bug_test2.Rdata")
	

	# Individual regeneration times
	regen_times = (which(regens == 1)) -
	         append(c(0), which(regens == 1)[-length(which(regens == 1))])

	mu = mean(regen_times)

	T = append(c(0),which(regens == 1))
	Z = matrix(0, nrow = sum(regens), ncol = 2)

	for(i in 1:(length(T)-1)){
		for(j in (T[i] + 1):T[i+1]){
			Z[i,1] = Z[i,1] + xydata[j,1]
			Z[i,2] = Z[i,2] + xydata[j,2]
		}
	}

	Z_bar = colMeans(Z)

	W = Z - (regen_times %*% t(Z_bar))/mu
	#ess = multiESS(xydata[,-3], covmat = variance(Z, regen_times))
	#regVAR1 = norm(Tr - variance(Z, regen_times), type = "F")
	cov = variance(W, regen_times)
	return(list("cov" = cov))
}

Gen_data = function(nsim){
	x = rep(0, nsim)
	y = rep(0, nsim)
	s = rep(0, nsim)
	x[1] = 1 # 1st initial point
	y[1] = rnorm(1, mu2 + rho*(x[1] - mu1), sqrt(omega2 - rho^(2)/omega1))
	if(y[1] - y_star > 0){
			s[1] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
		    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
		    exp(- rho^(2)*((y[1] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
		    exp(-h*rho*omega2*(y[1] - y_star)/(omega1*omega2 - rho^(2)))
		}else{
			s[1] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
		    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
		    exp(- rho^(2)*((y[1] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
		    exp(h*rho*omega2*(y[1] - y_star)/(omega1*omega2 - rho^(2)))
		}

	for(i in 2:nsim){
		x[i] = rnorm(1, mu1 + rho*(y[i-1] - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))
		y[i] = rnorm(1, mu2 + rho*(x[i] - mu1)/omega1, sqrt(omega2 - rho^(2)/omega1))
		if(y[i] - y_star > 0){
			s[i] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
		    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
		    exp(- rho^(2)*((y[i] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
		    exp(-h*rho*omega2*(y[i] - y_star)/(omega1*omega2 - rho^(2)))
		}else{
			s[i] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
		    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
		    exp(- rho^(2)*((y[i] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
		    exp(h*rho*omega2*(y[i] - y_star)/(omega1*omega2 - rho^(2)))
		}	
	}
	return(cbind(x, y, s))
}

