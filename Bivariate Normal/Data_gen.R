set.seed(1234)

load("Expected_S.Rdata")
Es
library(mcmcse)

# Initialization
nsim = 1e7

omega1 = 1
omega2 = 5
mu1 = 0
mu2 = 0
h = 1
y_star = 0
rho = 0.05

log_10_n = c(4,5,6,7)


Gen_data = function(nsim){
	x = rep(0, nsim)
	y = rep(0, nsim)
	s = rep(0, nsim)
	x[1] = 5 # 1st initial point
	y[1] = rnorm(1, mu2 + rho*(x[1] - mu1), sqrt(omega2 - rho^(2)/omega1))
	s[1] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
		    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
	        exp(- rho^(2)*((y[1] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
	         exp(-h*rho*omega2*(y[1] - y_star)/(omega1*omega2 - rho^(2)))


	for(i in 2:nsim){
		x[i] = rnorm(1, mu1 + rho*(y[i-1] - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))
		y[i] = rnorm(1, mu2 + rho*(x[i] - mu1), sqrt(omega2 - rho^(2)/omega1))
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
	return(list("x" = x, "y" = y, "s" = s))
}




# True Asymptotic Variance
Tr = matrix(c(omega1*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2)), 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)), 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)), omega2*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2))), nrow = 2)
Tr
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