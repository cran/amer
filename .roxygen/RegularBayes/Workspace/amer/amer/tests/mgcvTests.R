# Compare fits and fitting times for amer and mgcv
# data generating code taken from example(gamm) 
#
# Author: fabians
###############################################################################
if(FALSE)
{

detach(package:lme4)
require(nlme)
detach(nlme)
require(lme4)



library(mgcv)

## AMM
dat <- gamSim(6,n=4000,scale=.2)
detach(package:lme4)
library(nlme)
system.time(b2<-gamm(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat,random=list(fac=~1)))
plot(b2$gam,pages=1)
summary(b2$lme)
mean((fitted(b2$lme, level=3) - dat$f)^2)


detach(package:nlme)
library(lme4)
system.time(m2 <- amer(y ~ tp(x0, degree=2, k=10) + tp(x1, degree=2, k=10) + tp(x2, degree=2, k=10) + tp(x3, degree=2, k=10) + (1|fac), data =dat))
summary(m2)
mean((fitted(m2)-dat$f)^2)
getOption("device")() 
plotF(m2, addConst=F)

getOption("device")() 
par(mfrow=c(2,2))
with(dat,{
		plot(x0[o], (f0)[o])
		plot(x1[o], (f1)[o])
		plot(x2[o], (f2)[o])
		plot(x3[o], (f3)[o])
	}
)

plotF(m2, which="f.x0", addConst=F)
o <- order(dat$x0)
with(dat, lines(x0[o], (f0)[o]))


dat <- gamSim(6,n=400,scale=.2,dist="poisson")
detach(package:lme4)
library(nlme)
system.time(b3<-gamm(y~s(x0)+s(x1)+s(x2)+s(x3),family=poisson,
		data=dat,random=list(fac=~1)))
plot(b3$gam,pages=1)
summary(b3$lme)
mean((fitted(b3$lme, level=0) - log(dat$f))^2)


detach(package:nlme)
library(lme4)
system.time(m3 <- amer(y ~ tp(x0, degree=1, k=10) + tp(x1, degree=1, k=10) + tp(x2, degree=1, k=10) + tp(x3, degree=1, k=10) + (1|fac), family=poisson, data =dat))
summary(m3)
mean((log(fitted(m3))-log(dat$f))^2)
getOption("device")() 
plotF(m3, addConst=F)

}
