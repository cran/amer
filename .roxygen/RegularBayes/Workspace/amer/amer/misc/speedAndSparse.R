# compare speeds or different sparsity of base
#
# Author: fabians
###############################################################################

		
library(lme4)
library(mgcv)
sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(trace) cat(nm,":")           
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}
sourceDir("F:\\lme4Spline\\amer\\R")

set.seed(1234)

n <- 5000
snr <- .5
x <- sort(runif(n))
f <- dbeta(x,2,5)
y <-  f + sqrt(var(f)/snr)*rnorm(n)
d <- data.frame(y, x)


reps <- 10

## B-splines with diff. order = 0 ( -> Z is a band matrix)
fivenum(replicate(10, system.time(bs0 <<- amer(y ~ bsp(x, k=25, diff.ord = 0), data = d))[3]))
# 0.86    0.87    0.88    0.93    1.08 
mean(bs0@Zt==0)
#0.84

## constant TP
fivenum(replicate(10, system.time(tp0 <<- amer(y ~ tp(x, k=25, degree = 0), data = d))[3]))
#   1.78    1.81    1.81    1.81    2.02 
mean(tp0@Zt==0)
#0.52


## quadratic TP
fivenum(replicate(10, system.time(tp2 <<- amer(y ~ tp(x, k=25, degree = 2), data = d))[3]))
# 1.540   1.570   1.595   1.620   2.070 
mean(tp2@Zt==0)
#0.52

## quadratic TP with penalized global quadratic
fivenum(replicate(10, system.time(tpU2 <<- amer(y ~ tpU(x, k=25), data = d))[3]))
#1.610   1.630   1.635   1.660   1.990 
mean(tpU2@Zt==0)
#0.498

## B-splines with diff. order = 2 ( -> Z has no zeroes)
fivenum(replicate(10, system.time(bs3 <<- amer(y ~ bsp(x, k=25), data = d))[3]))
#3.810   3.920   4.025   4.220   5.300 
mean(bs3@Zt==0)
#0


######
# - no speed penalty for additional fixef (see tp0 vs tp2)
# - sparsity matters!

graphics.off()
plotF(bs0); lines(x, f, col=2)
windows(); plotF(tp0); lines(x, f, col=2)
windows(); plotF(tp2); lines(x, f, col=2)
windows(); plotF(tpU2); lines(x, f, col=2)
windows(); plotF(bs3); lines(x, f, col=2)

library(gamm4)
fivenum(replicate(10, system.time(cr4 <<- gamm4(y ~ s(x, k=25, bs="cr"), data = d))[3]))
#  3.030   3.170   3.395   3.820   5.030 
windows(); plot(cr4$gam); lines(x, f-mean(f), col=2)

detach(package:lme4)
fivenum(replicate(10, system.time(cr <<- gamm(y ~ s(x, k=25, bs="cr"), data = d))[3]))
# 29.72   30.78   31.30   32.01   33.40 
windows(); plot(cr$gam); lines(x, f-mean(f), col=2)
