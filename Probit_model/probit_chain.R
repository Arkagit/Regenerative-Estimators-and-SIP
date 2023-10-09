library(MASS)
library(truncnorm)


## AC for Bayesian probit regression
probit_gibbs <- function(dat, nsim){
  # Separate covariates and response
  x <- as.matrix(dat[, - dim(dat)[2]])
  p <- dim(x)[2]
  y <- dat[, dim(dat)[2]]  
  m <- length(y)
  
  #return objects  
  betas <- matrix(0, nrow = nsim, ncol = p)  
  z <- matrix(0, nrow = nsim, ncol = m)
  
  # starting value of beta are the mle estaimtes  
  z.curr <- numeric(length = m)  
  fit <- glm(y ~ -1 + x, family = binomial(link = "probit") )  
  beta.curr <- coef(fit)

  # Needed for the truncated normal  
  lower <- numeric(length = m)  
  upper <- numeric(length = m)

  
  for(i in 1:m)
  {
    if(y[i] == 0)
    {
      lower[i] <- -Inf      
      upper[i] <- 0
      
    }else
    {
      lower[i] <- 0      
      upper[i] <- Inf
    }
  }
  z.curr <- rtruncnorm(m, a = lower, b = upper, mean = x%*%beta.curr, sd = 1)
  
  errors <-  mvrnorm(n = nsim, rep(0,p), inv.xtx)
  
  for(i in 1:nsim){
    #if(i %% (nsim/10) == 0) print(paste(i/nsim * 100, "% done"))
    
    comp <- rbinom(1, size = 1, prob = 1/2)
    if(comp == 1) z.curr <- rtruncnorm(m, a = lower, b = upper, mean = x%*%beta.curr, sd = 1)
    if(comp == 0) beta.curr <- inv.xtx.tx %*% z.curr + errors[i, ]
    
    betas[i, ] <- beta.curr
    z[i, ] <- z.curr
  }

  return(list("beta" = betas, "z" = z))
}




###############################################################################
###############################################################################
# Required Functions
# Density for MV Normal
density_mv = function(input, mu, V){
  lngt = length(input)
  obj = exp(-0.5*t(input - mu)%*%solve(V)%*%as.matrix(input - mu))/sqrt(det(V)*((2*pi)^(lngt)))
  return(obj)
}

# Density for Truncated Normal
density_tn = function(inp, betaj){
  obj1 = numeric(length = m)
  for (t in 1:m) {
    if(inp[t] > 0){
      obj1[t] = dtruncnorm(inp[t], a = 0, b = Inf, mean = (as.matrix(x)%*%betaj)[t], 1)
    }else{
      obj1[t] = dtruncnorm(inp[t], a = -Inf, b = 0, mean = (as.matrix(x)%*%betaj)[t], 1)
    }
  }
  return(obj1)
}




# Function of success probabilities for Cuoidal D* for Random Scan
etaj_cube_randm2 <- function(betaj0, zj0, betaj2, zj2){
  ret <- 0
  t.starj <- t(zj0 - z.star)%*%as.matrix(x)
  
  if(prod((betaj1 < d & betaj1 > c)) == 1){

    ret <- sum(t.starj * (c* (t.starj > 0) + d * (t.starj < 0)))- 
            0.5*t(zj0)%*%as.matrix(x)%*%inv.xtx%*%t(x)%*%zj0 +
            0.5*t(z.star)%*%as.matrix(x)%*%inv.xtx%*%t(x)%*%z.star + 
            log(density_mv(betaj1, inv.xtx.tx %*% z.star, inv.xtx)) + sum(log(density_tn(zj1, betaj1)))
    
    ret <- ret - log(density_mv(betaj1, inv.xtx.tx%*%zj0, inv.xtx) + 
          prod(density_tn(zj1,betaj0)) + 
          density_mv(betaj2, inv.xtx.tx%*%zj2, inv.xtx)*prod(density_tn(zj2,betaj0)) +
            density_mv(betaj2, inv.xtx.tx%*%zj0, inv.xtx)*prod(density_tn(zj2,betaj2)) )
    
    ret <- exp(ret)
  }
  return(ret)
}

# Function of success probabilities for Cuoidal D* for Deterministic Scan
etaj_cube_det <- function(betaj0, zj0, betaj1, zj1){
  ret <- 0
  t.starj <- t(zj0 - z.star)%*%as.matrix(x)
  if(prod((betaj1 < d & betaj1 > c)) == 1){
    ret <- sum(t.starj * (c* (t.starj > 0) + d * (t.starj < 0) - betaj1))
    ret <- exp(ret)
  }
  return(ret)
}

# Regenerative Variance
#regen_var <- function(Z, mu, regen_steps){
#  p = dim(Z)[2]
#
#  A1 = matrix(0, nrow= p, ncol = p)
#  A2 = matrix(0, nrow= p, ncol = p)
#  A3 = matrix(0, nrow= p, ncol = p)
#
#  A1 = var(Z)
# for (i in 1:(length(regen_steps)-1)) {
 #   A2 = A2 + (Z[i,] - Z_bar)%*%t(Z[i+1,] - Z_bar)/dim(Z)[1]
  #  A3 = A3 + (Z[i+1,] - Z_bar)%*%t(Z[i,] - Z_bar)/dim(Z)[1]
 # }

#  Sigma_1 = (A1 + A2 + A3)/mu

 # return(Sigma_1)
#}

regen_var <- function(Z, mu, regen_steps){
  Z_bar = colMeans(Z)
  p = dim(Z)[2]
  betak = Z_bar/mu

  regen_steps = which(regen_det == 1) + 1
  regen_points = c(1, regen_steps[- length(regen_steps)])
  regen_length = regen_steps - regen_points

  # Calculating sample variance using regenerations
  W = matrix(0, nrow = length(regen_steps), ncol = p)

  for (j in 1:length(regen_steps)) {
    W[j, ] = (Z[j, ] - t(Z_bar)*regen_length[j]/tau_bar)
  }

  cov_lag1 = matrix(0, nrow = p, ncol = p) # Initializing auto-covariances
  for (i in 1:(length(regen_steps) - 1)) {
    cov_lag1 = cov_lag1 + W[i, ]%*%t(W[i + 1, ])/length(regen_steps)
  }
  Sigma_2 = (var(W) + cov_lag1 + t(cov_lag1)) / mu
  return(Sigma_2)
}

regen_var2 <- function(Z, mu, regen_steps){
  Z_bar = colMeans(Z)
  p = dim(Z)[2]

  regen_steps = which(regen_det == 1) + 1
  regen_points = c(1, regen_steps[- length(regen_steps)])
  regen_length = regen_steps - regen_points

  A = matrix(0, nrow = p, ncol = p)

  for (j in 1:length(regen_steps)) {
    A = A + (Z[j,] - beta_reg_mean*regen_length[j])%*%t(Z[j,] - beta_reg_mean*regen_length[j])
  }
  return(A/(length(regen_steps)*mu))
}



regen_var3 <- function(Z, mu, regen_steps){
  Z_bar = colMeans(Z)
  p = dim(Z)[2]
  betak = Z_bar/mu

  regen_steps = which(regen_det == 1) + 1
  regen_points = c(1, regen_steps[- length(regen_steps)])
  regen_length = regen_steps - regen_points
  
  # Calculating sample variance using regenerations
  W = matrix(0, nrow = length(regen_steps), ncol = p)

  for (j in 1:length(regen_steps)) {
    W[j, ] = (Z[j, ] - t(Z_bar))  - t(betak * (regen_length[j]  - mu))
  }

  Sigma_3 = mcse.multi(W, method = "bm", size = 2)$cov/(mu)
  return(Sigma_3)
}


