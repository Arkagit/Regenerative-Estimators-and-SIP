source("Initialization.R")
load("NE.Rdata")
     
log_10_n = log(samp_size)/log(10)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


#pdf("boxplot.pdf", height = 6, width = 6)
#par(mar = c(5.1, 4.8, 4.1, 2.1))
#boxplot(ess_reg1,ess_reg2, ess_bmopt, ess_bmth, xlab = "Variance Estimators", ylab = "Effective Sample Size",
#	main = "ESS for different variance estimators", col = c("blue", "green","yellow","purple"),
 #  names = c("1-step","2-step","BM(opt)","BM(theory)"))
#dev.off()

BMopt_frob = apply(norm_bmopt, 2, sd)/sqrt(reps)
BMth_frob = apply(norm_bmth, 2, sd)/sqrt(reps)
regVAR1_frob = apply(norm_reg1, 2, sd)/sqrt(reps)
regVAR2_frob = apply(norm_reg2, 2, sd)/sqrt(reps)


pdf("plot_frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(log_10_n, colMeans(norm_reg2), type = 'l', ylim = c(0, 6), xlab = expression(log[10]~n), ylab = "Frobenius norm",
	main = "Frobenius distance from True Covariance")
segments(x0 = log_10_n, y0 = colMeans(norm_reg2) - 1.96*regVAR2_frob, 
	y1 = colMeans(norm_reg2) + 1.96*regVAR2_frob)

lines(log_10_n, colMeans(norm_reg1), type = 'l', col = "red")
segments(x0 = log_10_n, y0 = colMeans(norm_reg1) - 1.96*regVAR1_frob, 
	y1 = colMeans(norm_reg1) + 1.96*regVAR1_frob, col = "red")

lines(log_10_n, colMeans(norm_bmopt), type = 'l', col = "blue")
segments(x0 = log_10_n, y0 = colMeans(norm_bmopt) - 1.96*BMopt_frob, 
	y1 = colMeans(norm_bmopt) + 1.96*BMopt_frob, col = "blue")

lines(log_10_n, colMeans(norm_bmth), type = 'l', col = "skyblue")
segments(x0 = log_10_n, y0 = colMeans(norm_bmth) - 1.96*BMth_frob, 
	y1 = colMeans(norm_bmth) + 1.96*BMth_frob, col = "skyblue")

abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("2-step", "1-step", "BM(optimal)", "BM(Theory)"), 
	col = c("black", "red","blue", "skyblue"), lty = 1, cex=0.5)

dev.off()


BMopt_ess = apply(ess_bmopt, 2, sd)/sqrt(reps)
BMth_ess = apply(ess_bmth, 2, sd)/sqrt(reps)
regVAR1_ess = apply(ess_reg1, 2, sd)/sqrt(reps)
regVAR2_ess = apply(ess_reg2, 2, sd)/sqrt(reps)
ymax = max(ess_reg1, ess_reg2, ess_bmth, ess_bmopt)

pdf("plot_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(log_10_n, colMeans(ess_reg2), type = 'l', ylim = c(0.73, 0.77),xlab = expression(log[10]~n), ylab = "ESS/n",
	main = "Effective Sample Size")
segments(x0 = log_10_n, y0 = colMeans(ess_reg2) - 1.96*regVAR2_ess, 
	y1 = colMeans(ess_reg2) + 1.96*regVAR2_ess)

lines(log_10_n, colMeans(ess_reg1), type = 'l', col = "red")
segments(x0 = log_10_n, y0 = colMeans(ess_reg1) - 1.96*regVAR1_ess, 
	y1 = colMeans(ess_reg1) + 1.96*regVAR1_ess, col = "red")

lines(log_10_n, colMeans(ess_bmopt), type = 'l', col = "blue")
segments(x0 = log_10_n, y0 = colMeans(ess_bmopt) - 1.96*BMopt_ess, 
	y1 = colMeans(ess_bmopt) + 1.96*BMopt_ess, col = "blue")

lines(log_10_n, colMeans(ess_bmth), type = 'l', col = "skyblue")
segments(x0 = log_10_n, y0 = colMeans(ess_bmth) - 1.96*BMth_ess, 
	y1 = colMeans(ess_bmth) + 1.96*BMth_ess, col = "skyblue")

abline(h = (det(Target_mat)/det(Tr))^(1/2), lty = 2)
legend("topright", bty = "n",legend = c("2-step", "1-step", "BM(cuberoot)", "BM(sqroot)", "TRUE"), 
	col = c("black", "red","blue", "skyblue"), lty = c(1,1,1,1,2), cex=0.5)

dev.off()
