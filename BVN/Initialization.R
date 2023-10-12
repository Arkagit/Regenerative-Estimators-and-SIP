# Reps
nsim = 1e8
nsize = 1e7
reps = 500
samp_size = c(1e4, 5e4, 1e5, 3e5, 5e5, 8e5, 1e6)

# Study parameters
omega1 = 1
omega2 = 5
rho = 1.5
mu1 = 0
mu2 = 0
y_star = 0
h = 1
Target_mat = matrix(c(omega1, rho, rho, omega2), nrow = 2, ncol = 2)

# Output initialization
x = rep(0, nsim)
y = rep(0, nsim)
s = rep(0, nsim)


# True Asymptotic Variance
Tr = matrix(c(omega1*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2)),
 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)), 
 2*omega1*omega2*rho/(omega1*omega2 - rho^(2)),
  omega2*(omega1*omega2 + rho^(2))/(omega1*omega2 - rho^(2))), nrow = 2); Tr