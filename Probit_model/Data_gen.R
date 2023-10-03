## Entering Data

#x1.dummy <- seq(-3, 1.5, by = 0.5)
#x2.dummy <- seq(0, 2, by = 0.5)

## The two column covariates
#x1 <- rep(x1.dummy, each = length(x2.dummy))
#x2 <- rep(x2.dummy, times = length(x1.dummy))

# Total number of cases. Follow the rows in the image above

#n <- c(1, rep(0,4), 3, rep(0,4), 7, rep(0,3),
#       1, 6, 1, rep(0,3), 6,1,1,0,1,4,0,0,1,0,3,
#       0,1,1,0,4,0,1,1,1,1,0,1,1,4,1,0,0,2,0)

# Number of Lupus cases, follow the numerator across rows in the image

#y <- c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,3,1,1,1,1,1,1,4,1,2)

# Following code makes a data matrix of all cases and lupus cases removing
# the combinations of covariates with noncases

#dat.binom <- cbind(n,rep(1, 50), x1, x2)
#dat.binom <- dat.binom[!(dat.binom[,1] == 0), ]
#dat.binom <- cbind(dat.binom,y)

# Following code essentially converts binomial data to bernoulli data
#dummy <- sapply(1:length(y), function(x) matrix(rep(dat.binom[x, -1], dat.binom[x,1]),nrow = dat.binom[x,1], byrow = TRUE) )
#dat <- do.call(rbind, dummy)
# dat <- dat[55,]
#dat[(dat[,4] >0),4] <- 1

# Fixing a mistake

#dat[39,4] <- 0
foo <- read.csv("titanic.csv")
head(foo)
x = as.matrix(foo[,-1])
y = foo[,1]
dat = cbind(x,y)
#colnames(dat) <- c("intercept", "x1", "x2", "y")
dim(dat)
#x <- dat[, -4]

# Summary Statistics Needed Later
xtx <- t(x)%*%x
  
inv.xtx <- solve(xtx)
inv.xtx.tx <- solve(xtx) %*% t(x)
x <- dat[, -4]
p <- dim(x)[2]
y <- dat[, 4]  
m <- length(y)