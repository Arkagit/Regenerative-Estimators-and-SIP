nsim = 1e8
x = rep(0, nsim)
y = rep(0, nsim)

omega1 = 1
omega2 = 5
mu1 = 0
mu2 = 0
h = 1
y_star = 0
rho = 0.05
s = rep(0, nsim)

pr_D = pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
	    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))

p_xgy = function(ndat){
	x = rep(1000, ndat)
	y = rep(1000, ndat)
	for(i in 1:ndat){
		while(abs(x[i]) > h){
			x[i] = rnorm(1, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))
		}
		y[i] = rnorm(1, mu2 + rho*(x[i] - mu1), sqrt(omega2 - rho^(2)/omega1))
		if(i/1e5 %% 1 == 0){
			print(i)
		}
	}
	return(cbind(x, y))
}


pilot_data = p_xgy(nsim)
head(pilot_data)

for(i in 2:nsim){
	if(pilot_data[i,2] - y_star > 0){
		s[i] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
	    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
	    exp(- rho^(2)*((pilot_data[i,2] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
	    exp(-h*rho*omega2*(pilot_data[i,2] - y_star)/(omega1*omega2 - rho^(2)))
	}else{
		s[i] = (pnorm(h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2)) - 
	    pnorm(-h, mu1 + rho*(y_star - mu2)/omega2, sqrt(omega1 - rho^(2)/omega2))) * 
	    exp(- rho^(2)*((pilot_data[i,2] - mu2)^(2) - (y_star - mu2)^(2)) / (2*omega2*(omega1*omega2 - rho^(2)))) * 
	    exp(h*rho*omega2*(pilot_data[i,2] - y_star)/(omega1*omega2 - rho^(2)))
	}	
}

Es = mean(s); Es
save("Es" = Es, file = "Expected_S.Rdata")

