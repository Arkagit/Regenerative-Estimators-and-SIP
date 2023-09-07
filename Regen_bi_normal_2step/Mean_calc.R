nsim = 1e7
x = rep(0, nsim)
y = rep(0, nsim)
s = rep(0, nsim)

y_star = 0

pr_D = pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))

p_xgy = function(ndat){
	x = rep(1000, ndat)
	y = rep(1000, ndat)
	for(i in 1:ndat){
		while(abs(x[i]) > 1){
			x[i] = rnorm(1, 2*y_star/5, 1/sqrt(5))
		}
		y[i] = rnorm(1, 2*x[i], 1)
	}
	return(cbind(x, y))
}


pilot_data = p_xgy(nsim)

for(i in 2:nsim){
	if(pilot_data[i,2] > 0){
		s[i] = (pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))) * exp(-(2*pilot_data[i,2]^{2} + 50*pilot_data[i,2])/5)
	}else{
		s[i] = (pnorm(1, 0, 1/sqrt(5)) - pnorm(-1, 0, 1/sqrt(5))) * exp(-(2*pilot_data[i,2]^{2} - 50*pilot_data[i,2])/5)
	}	
}

Es = mean(s)
Es <- Es; Es

