source("Data_gen.R")
load("1step.Rdata")
load("2step.Rdata")

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}



BMopt_frob = apply(BMopt, 2, sd)/sqrt(reps)
BMsq_frob = apply(BMsq, 2, sd)/sqrt(reps)
regVAR1_frob = apply(regVAR1, 2, sd)/sqrt(reps)
regVAR2_frob = apply(regVAR2, 2, sd)/sqrt(reps)


pdf("plot_frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(log_10_n, colMeans(regVAR2), type = 'l', ylim =  c(-1, 25), xlab = "log_10_(n)", ylab = "Frobenius norm",
	main = "Frobenius distance from True Covariance")
segments(x0 = log_10_n, y0 = colMeans(regVAR2) - 1.96*regVAR2_frob, 
	y1 = colMeans(regVAR2) + 1.96*regVAR2_frob)

lines(log_10_n, colMeans(regVAR1), type = 'l', col = "red")
segments(x0 = log_10_n, y0 = colMeans(regVAR1) - 1.96*regVAR1_frob, 
	y1 = colMeans(regVAR1) + 1.96*regVAR1_frob, col = "red")

lines(log_10_n, colMeans(BMopt), type = 'l', col = "blue")
segments(x0 = log_10_n, y0 = colMeans(BMopt) - 1.96*BMopt_frob, 
	y1 = colMeans(BMopt) + 1.96*BMopt_frob, col = "blue")

lines(log_10_n, colMeans(BMsq), type = 'l', col = "skyblue")
segments(x0 = log_10_n, y0 = colMeans(BMsq) - 1.96*BMsq_frob, 
	y1 = colMeans(BMsq) + 1.96*BMsq_frob, col = "skyblue")

abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("2-step", "1-step", "BM(optimal)", "BM(Theory)"), 
	col = c("black", "red","blue", "skyblue"), lty = 1, cex=0.5)

dev.off()
