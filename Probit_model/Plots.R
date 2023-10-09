library(mcmcse)
load("variance.Rdata")
     
log_10_n = log(samp_size)/log(10)

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

norm_bmopt = matrix(0, nrow = reps, ncol = length(samp_size))
norm_bmth = matrix(0, nrow = reps, ncol = length(samp_size))

for(i in 1:length(reps)){
	for(j in 1:length(samp_size)){
		norm_bmth[i,j] = norm(EST[[i]][[j]][[1]], type = "F")
		norm_bmopt[i,j] = norm(EST[[i]][[j]][[2]], type = "F")
	}
}
#pdf("boxplot.pdf", height = 6, width = 6)
#par(mar = c(5.1, 4.8, 4.1, 2.1))
#boxplot(ess_reg1,ess_reg2, ess_bmopt, ess_bmth, xlab = "Variance Estimators", ylab = "Effective Sample Size",
#	main = "ESS for different variance estimators", col = c("blue", "green","yellow","purple"),
 #  names = c("1-step","2-step","BM(opt)","BM(theory)"))
#dev.off()

BMopt_frob = apply(norm_bmopt, 2, sd)/sqrt(reps)
BMth_frob = apply(norm_bmth, 2, sd)/sqrt(reps)


pdf("plot_frob.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(log_10_n, colMeans(norm_bmopt), ylim = c(-0.002, 0.003), type = 'l', col = "blue",
 ylab = "Frobenius Norm", main = "Frobenius Norm for Batch Means estimator")
segments(x0 = log_10_n, y0 = colMeans(norm_bmopt) - 1.96*BMopt_frob, 
	y1 = colMeans(norm_bmopt) + 1.96*BMopt_frob, col = "blue", lty = 2)

lines(log_10_n, colMeans(norm_bmth), type = 'l', col = "skyblue")
segments(x0 = log_10_n, y0 = colMeans(norm_bmth) - 1.96*BMth_frob, 
	y1 = colMeans(norm_bmth) + 1.96*BMth_frob, col = "skyblue", lty = 2)

#abline(h = 0, lty = 2)
legend("topright", bty = "n",legend = c("BM(optimal)", "BM(Theory)"), 
	col = c("blue", "skyblue"), lty = 2, cex=0.5)

dev.off()

########################################################################
ess_bmopt = matrix(0, nrow = reps, ncol = length(samp_size))
ess_bmth = matrix(0, nrow = reps, ncol = length(samp_size))

for(i in 1:length(reps)){
	for(j in 1:length(samp_size)){
		ess_bmth[i,j] = EST[[i]][[j]][[3]]/samp_size[j]
		ess_bmopt[i,j] = EST[[i]][[j]][[4]]/samp_size[j]
	}
}


BMopt_ess = apply(ess_bmopt, 2, sd)/sqrt(reps)
BMth_ess = apply(ess_bmth, 2, sd)/sqrt(reps)

pdf("plot_ess.pdf", height = 6, width = 6)
par(mar = c(5.1, 4.8, 4.1, 2.1))

plot(log_10_n, colMeans(ess_bmopt), type = 'l',  col = "blue", ylim = c(-0.0002, 0.0005),
	ylab = "ESS/n", main = "Effective Sample size of Batch Means estimator")
segments(x0 = log_10_n, y0 = colMeans(ess_bmopt) - 1.96*BMopt_ess, 
	y1 = colMeans(ess_bmopt) + 1.96*BMopt_ess, col = "blue", lty = 2)

lines(log_10_n, colMeans(ess_bmth), type = 'l', col = "skyblue")
segments(x0 = log_10_n, y0 = colMeans(ess_bmth) - 1.96*BMth_ess, 
	y1 = colMeans(ess_bmth) + 1.96*BMth_ess, col = "skyblue", lty = 2)

#abline(h = (det(Target_mat)/det(Tr))^(1/2), lty = 2)
legend("topright", bty = "n",legend = c("BM(cuberoot)", "BM(sqroot)"), 
	col = c("blue", "skyblue"), lty = 2, cex=0.5)

dev.off()