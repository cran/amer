#2D-tensor spline basis with marginal TP-basis functions
# 
# Author: fabians
###############################################################################
source("F:\\lme4Spline\\amer\\R\\reparameterizeDesign.R")
source("F:\\lme4Spline\\amer\\R\\tpU2d.R")
library(amer)
library(mgcv)
library(mvtnorm)

sqrt.n <- 50
n <- sqrt.n^2
snr <- 100

grid <- expand.grid(seq(-1, 1, l=sqrt.n), 
		seq(-1, 1, l=sqrt.n))
x <- grid[,1]
z <- grid[,2]


rho <-.95
m <- c(0,0)
f <- dmvnorm(grid, mean=m, sigma= matrix(c(1, rho, rho, 1),2,2))+dmvnorm(grid, mean=-m, sigma= matrix(c(1, -rho, -rho, 1),2,2))
#rotsin <- function(m) {
#	x<-m[,1]
#	y<-m[,2]
#	r <- sqrt(x^2+y^2)
#	return(10 * sin(r)/r)
#}
#f <- rotsin(10*grid)


y <-  f + sqrt(var(f)/snr)*rnorm(n)
persp(matrix(f,sqrt.n,sqrt.n))

d <- data.frame(y,x,z)

system.time(test1 <- amer(y~ tpU2d(x,z), data=d, basisGenerators=c("tpU2d")))
getF(test1)


system.time(ref <- gam(y~ te(x,z, k=c(7,7)), data=d))

layout(t(1:4))
persp(matrix(y,sqrt.n,sqrt.n))
persp(matrix(f,sqrt.n,sqrt.n))
persp(matrix(fitted(test1),sqrt.n,sqrt.n))
persp(matrix(fitted(ref),sqrt.n,sqrt.n))



##################################
	## sanity check: 
	x11(50,50)
	par(mar= par()$mar/4)
	layout(matrix(1:((kx+1)*(kz+1)), kx+1,kz+1))
	plot(0)
	for(i in 2:((kx+1)*(kz+1))) persp(sort(unique(x)), sort(unique(z)), matrix(D[,i],sqrt.n,sqrt.n,byrow=F), main=paste(i,":",diag(P)[i]))







