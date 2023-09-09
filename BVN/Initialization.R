# Reps
nsim = 1e8
nsize = 1e7
reps = 10
samp_size = c(1e4, 1e5, 1e6, 1e7)

# Study parameters
omega1 = 1
omega2 = 5
rho = 0.05
mu1 = 0
mu2 = 0
y_star = 0
h = 1

# Output initialization
x = rep(0, nsim)
y = rep(0, nsim)
s = rep(0, nsim)


# True Asymptotic Variance
Tr = matrix(c(omega1*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2)),
 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)), 
 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)),
  omega2*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2))), nrow = 2); Tr