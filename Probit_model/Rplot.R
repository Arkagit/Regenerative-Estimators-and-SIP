set.seed(4321)

library(MASS)

library(Matrix)

library(truncnorm)


source("Small_set.R")



pdf(file="Rplot.pdf")
#ACF of beta samples from AC
par(mfrow = c(2,3))
acf(out.1e6$beta[,1], ylim = c(0,1), main = expression(paste("AC for ", beta(1))))
acf(out.1e6$beta[,2], ylim = c(0,1), main = expression(paste("AC for ", beta(2))))
acf(out.1e6$beta[,3], ylim = c(0,1), main = expression(paste("AC for ", beta(3))))
acf(out.1e6$beta[,4], ylim = c(0,1), main = expression(paste("AC for ", beta(4))))
acf(out.1e6$beta[,5], ylim = c(0,1), main = expression(paste("AC for ", beta(5))))
acf(out.1e6$beta[,6], ylim = c(0,1), main = expression(paste("AC for ", beta(6))))
dev.off()