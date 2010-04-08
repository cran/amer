# Basic tests of amer-functionality:
#
# Author: fabians
###############################################################################
library(amer)

#pdf("optionsTests.pdf")


###############################
# simple additive model
###############################
set.seed(1234)
n <- 50
snr <- 5
x <- sort(rexp(n))
x <- x/max(x)
f <- dbeta(x, 2, 5)
y <- f + sqrt(var(f)/snr) * rnorm(n)
d <- data.frame(y, x)

test00 <- amer(y ~ tp(x, degree = 0), data = d)
f00 <- plotF(test00, interval = "RW", type = "s")
lines(x, f, col = "grey")

test01 <- amer(y ~ tp(x, degree = 1), data = d)
f01 <- plotF(test01, n = 200, interval = "RW")
lines(x, f, col = "grey")

f01.MCMC <- getF(test01, n = 50, interval = "MCMC")
print(xyplot(attr(f01.MCMC, "mcmc")))

#sanity checks for MCMC:
f01.MCMCData <- as.data.frame(attr(f01.MCMC, "mcmc"))
unlist(c(fixef(test01), test01@ST, lme4:::sigma(test01)))
apply(f01.MCMCData, 2, quantile, probs = c(0.1, 0.25, 
				0.5, 0.75, 0.9), na.rm = T)
(bNormsObs <- unlist(lapply(ranef(test01), function(x) sqrt(sum(x^2)))))
bNormsMCMC <- apply(attr(f01.MCMC, "mcmc")@ranef, 
		2, function(x) return(sqrt(sum(x[1:14]^2, na.rm = T))))
quantile(bNormsMCMC, probs = c(0.1, 0.25, 0.5, 0.75, 
				0.9), na.rm = T)

test02 <- amer(y ~ tp(x, degree = 2), data = d)
f02 <- plotF(test02)
lines(x, f, col = "grey")

###############################
# by - option
###############################
set.seed(123554)
n <- 300
snr <- 5
x <- rbeta(n, 0.5, 1)
f1 <- dbeta(x, 2, 3)
f2 <- dbeta(x, 2, 5)
g1 <- factor(rep(letters[1:2], length = n))
f <- f1 * (g1 == "a") + (2 + f2) * (g1 == "b")
y <- f + sqrt(var(f)/snr) * rnorm(n)
z <- rnorm(n)
d <- data.frame(y, x, z, g1)

(test1 <- amer(y ~ g1 * z + tp(x, by = g1, k = 25, 
							degree = 2), data = d))
f1 <- plotF(test1, interval = "RW")
points(d$x, f, col = g1, pch = 19)
f11 <- getF(test1, newdata = d[order(d$x), ], interval = "MCMC")
lines(f11[[1]][[1]]$x, f11[[1]][[1]]$fhat, col = 3)
lines(f11[[1]][[1]]$x, f11[[1]][[1]]$lo, col = 3, 
		lty = 3)
lines(f11[[1]][[1]]$x, f11[[1]][[1]]$hi, col = 3, 
		lty = 3)
lines(f11[[1]][[2]]$x, f11[[1]][[2]]$fhat, col = 4)
lines(f11[[1]][[2]]$x, f11[[1]][[2]]$lo, col = 4, 
		lty = 3)
lines(f11[[1]][[2]]$x, f11[[1]][[2]]$hi, col = 4, 
		lty = 3)

(test12 <- amer(y ~ tp(x, by = g1, k = 25, degree = 2) + 
							tp(z, by = g1, k = 10, degree = 1), data = d))
f12 <- plotF(test12, addConst=c(TRUE, FALSE))


##########################
# allP - option
##########################
set.seed(12333554)
n <- 500
snr <- 100
grps <- 5
x <- rep(seq(0, 1, l = n/grps), times = grps)
g1 <- gl(n = grps, k = n/grps, labels = 1:grps)
i <- 2
f <- unlist(tapply(x, g1, function(x) {
					ff <- sqrt(i) + dbeta(x, 5, 2) + dbeta(x, 2 * sqrt(i), i/2)
					i <<- i + 1
					return(ff)
				}))
plot(x, f, col = g1)
y <- f + sqrt(var(f)/snr) * rnorm(n)
d <- data.frame(y, x = scale(x), g1)

#all penalized smooths
(test2 <- amer(y ~ tp(x, by = g1, k = 8, allPen = T), 
					data = d))
f2 <- plotF(test2, n = n/grps, interval = "none")
points(d$x, f, col = g1, pch = 19)
#points(d$x, fitted(test2), col=g1)

#check getF for supplied newdata
par(mfrow = c(1, 2))
f22 <- plotF(test2, n = n/grps, interval = "none", addConst = F, 
		auto.layout = F)
f2d <- getF(test2, newdata = d, addConst = F)
plot(0, 0, ylim = range(f22), xlim = range(d$x))
i <- 1
lapply(f2d[[1]], function(x) {
			points(x$x, x$fhat, col = i)
			i <<- i + 1
		})

#fit one global curve and penalized deviations for groups
(test21 <- amer(y ~ tp(x, k = 10, degree = 2) + tp(x, 
							by = g1, k = 15, degree = 1, allPen = T), data = d))
f21 <- plotF(test21, n = n/grps, int = "none", addConst = F)
fmat <- f21[[1]][[1]]$fhat + do.call(cbind, lapply(f21[[2]], 
				function(x) x$fhat)) + fixef(test21)["(Intercept)"]
plot(x, f, col = g1)
matlines(x[1:(n/grps)], fmat, type = "l", col = 1:5)

#fit with diagonal random effects
(test22 <- amer(y ~ tp(x, by = g1, k = 8, allPen = T, 
							diag = T), data = d))
f22 <- plotF(test22, n = n/grps, int = "none", lwd = 3)
plot(d$x, f, col = g1, pch = 19)
points(d$x, fitted(test22), col = g1)


########################
# varying - option
########################
set.seed(55423234)
n <- 100
snr <- 10
x <- sort(rbeta(n, 1, 2))
f1 <- dbeta(x, 5, 2)
z <- runif(n, -1, 1)
f <- f1 * z
y <- f + sqrt(var(f)/snr) * rnorm(n)
d <- data.frame(y, x, z)

(test3 <- amer(y ~ tp(x, k = 5, varying = z), data = d))
f3 <- plotF(test3, int = "RW")
f3 <- plotF(test3, int = "MCMC")
lines(x, f1, col = 2)
f32 <- getF(test3, newdata = d, int = "RW")
range(f - fitted(test3))

#this shouldn't work:
try(test31 <- amer(y ~ z + tp(x, varying = z), data = d))
try(test32 <- amer(y ~ poly(z, 2) + tp(x, varying = z), 
				data = d))

###############################
# generalized additive model
###############################
set.seed(12434)
n <- 400
x <- sort(runif(n, -1, 1))
f <- 2 * (scale(dnorm(x, m = -1, sd = 0.5)) + scale(dnorm(x, 
							m = 0.5, sd = 0.5)))
y <- rpois(n, exp(f))
d <- data.frame(y, x, f)
plot(x, jitter(y))
lines(x, exp(f))

(test4 <- amer(y ~ tp(x, degree = 2), data = d, family = poisson))
par(mfrow=c(1,2))
f4 <- plotF(test4, int = "RW", addConst = T, trans = exp)
lines(x, exp(f), col = "grey")
f42 <- plotF(test4, int = "RW", addConst = T)
lines(x, f, col = "grey") 

############################
# passing formula-objects:
############################
f <- as.formula("y ~ tp(x, degree = 2)")
(testForm <- amer(f, data = d, family = poisson))

#dev.off()

