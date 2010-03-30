if(FALSE){
options(digits=3, chmhelp=F, htmlhelp=T, width=120)
require(lme4)
sourceDir <- function(path, trace = TRUE, ...) {
	for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(trace) cat(nm,":")           
		source(file.path(path, nm), ...)
		if(trace) cat("\n")
	}
}
sourceDir("Z:\\lme4Spline\\amer\\R")






##############################
# predict function
##############################
n <- 600
grps <- 20
g <- gl(grps, n/grps)
bg <- bg0 <-  rnorm(grps)
f <- factor(sample(gl(2, n/2), rep=T))

x1 <- runif(n)
x2 <- runif(n)

#AM
truth <-  dbeta(x2, 4, 9)
y <- truth + .01*rnorm(n)
d1 <- data.frame(y, truth, x1, x2, f, g)

#AMM
truth <-  dbeta(x2, 4, 9) + rep(bg, e=n/grps)
y <- truth + .01*rnorm(n)
d2 <- data.frame(y, truth, x1, x2, f, g)

truth <-  dbeta(x2, 4, 9) + rep(bg, e=n/grps) + sqrt(abs(rep(bg, e=n/grps)))*x2
y <- truth + .01*rnorm(n)
d3 <- data.frame(y, truth, x1, x2, f, g)

truth <-  dbeta(x2, 4, 9) + rep(bg, e=n/grps) + sqrt(abs(rep(bg, e=n/grps)))*x2
y <- truth + rnorm(n)
d4 <- data.frame(y, truth, x1, x2, f, g)



newdata <- {
	g <- gl(grps, n/grps)
	bg <- rep(bg0, length=grps)
	f <- factor(sample(gl(2, n/2), rep=T))
	x1 <- runif(n)
	x2 <- runif(n)
	
	y1 <- truth1 <- dbeta(x2, 4, 9)
	y2 <- truth2 <-  dbeta(x2, 4, 9) + rep(bg, e=n/grps)
	y3 <- truth3 <-  dbeta(x2, 4, 9) + rep(bg, e=n/grps) + sqrt(abs(rep(bg, e=n/grps)))*x2
	y4 <- truth3 + rnorm(n)
			
	data.frame(y1, truth1, y2, truth2, y3, y4, truth3, x1, x2, f, g)
}	


# AM
(test1 <- amer(y ~ tp(x2), data = d1))
pred1 <- predict(test1, newdata)
plot(pred1$truth1, pred1$fit)
ggplot(aes(x=x2, y=y1), data=pred1) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))

#AMM
(test2 <- amer(y ~ tp(x2) + (1|g), data = d2))
summary(fitted(test2) - predict(test2, d2)$fit)
pred2 <- predict(test2, newdata)
plot(pred2$truth2, pred2$fit)
ggplot(aes(x=x2, y=y2), data=pred2) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))

(test3 <- amer(y ~ tp(x2) + (x2|g), data = d3))
pred3 <- predict(test3, newdata)
plot(pred3$truth3, pred3$fit)
ggplot(aes(x=x2, y=y3), data=pred3) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))


#AMM with linear effects
(test4 <- amer(y ~ x1 + f  + tp(x2) + (x2|g), data = d4))
pred4 <- predict(test4, newdata)
plot(pred4$truth3, pred4$fit)
ggplot(aes(x=x2, y=y3), data=pred3) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))


##########################
## additional grp-levels:
newdata <- {
	moreGrps <- 1.5
	g <- gl(moreGrps*grps, n/(moreGrps*grps))
	bg <- rep(bg0, length=moreGrps*grps)
	f <- factor(sample(gl(2, n/2), rep=T))
	x1 <- runif(n)
	x2 <- runif(n)
	
	y1 <- truth1 <- dbeta(x2, 4, 9)
	y2 <- truth2 <-  dbeta(x2, 4, 9) + rep(bg, e=n/(moreGrps*grps))
	y3 <- truth3 <-  dbeta(x2, 4, 9) + rep(bg, e=n/(moreGrps*grps)) + sqrt(abs(rep(bg, e=n/(moreGrps*grps))))*x2
	y4 <- truth3 + rnorm(n)
	
	data.frame(y1, truth1, y2, truth2, y3, y4, truth3, x1, x2, f, g)
}	
predMore <- predict(test2, newdata)
ggplot(aes(x=x2, y=y2), data=predMore) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))

## missing grp-levels:
newdata <- {
	lessGrps <- .5
	g <- gl(lessGrps*grps, n/(lessGrps*grps))
	bg <- bg0[1:(lessGrps*grps)]
	f <- factor(sample(gl(2, n/2), rep=T))
	x1 <- runif(n)
	x2 <- runif(n)
	
	y1 <- truth1 <- dbeta(x2, 4, 9)
	y2 <- truth2 <-  dbeta(x2, 4, 9) + rep(bg, e=n/(lessGrps*grps))
	y3 <- truth3 <-  dbeta(x2, 4, 9) + rep(bg, e=n/(lessGrps*grps)) + sqrt(abs(rep(bg, e=n/(lessGrps*grps))))*x2
	y4 <- truth3 + rnorm(n)
	
	data.frame(y1, truth1, y2, truth2, y3, y4, truth3, x1, x2, f, g)
}	
predLess <- predict(test2, newdata)
ggplot(aes(x=x2, y=y2), data=predLess) + geom_point(aes(group=g, col=g)) + facet_wrap(~g)+
		geom_line(aes(x=x2, y=fit, group=g, col=g))
}