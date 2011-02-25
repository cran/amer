#' Generate a reparameterized b-spline basis
#'
#' Generate a b-spline basis with equidistant knots in mixed model reparameterization
#' 
#' @param x covariate for the smooth function
#' @param k integer: dimensionality of the basis 
#' @param spline.degree integer: degree of B-splines (defaults to cubic)
#' @param diff.ord integer: order of the difference penalty on the un-reparamerized spline coefficients. Defaults to 2, that is, penalized deviations from linearity.
#' @param knots double: vector of (inner & outer) knot locations 
#' @param by factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.
#' @param allPen boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
#'                         if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X
#' @param varying numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying
#' @param diag logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}? 
#' @return list with entries:
#'			\code{"X"}: \code{n x (diff.ord - 1)} design matrix for unpenalized part (without intercept), 
#' 			\code{"Z"}: \code{n x (k - diff.ord+1)} design matrix for penalized part,
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{tp}}
bsp <- function(x, k=15, spline.degree = 3, diff.ord = 2, 
		knots=NULL, by=NULL, allPen = FALSE, varying = NULL, diag=FALSE)
{
	require(splines)
	call <- as.list(expand.call())
	stopifnot(diff.ord>=0, spline.degree>=0, is.numeric(x), k > spline.degree)
	
	if(is.null(call$knots)){
		knots.no <- k - spline.degree + 1 
		#generate a B-Spline-Matrix with equidistant knots (adapted from code by Thomas Kneib):
		n<-length(x)
		xl<-min(x)
		xr<-max(x)
		xmin<-xl-(xr-xl)/100
		xmax<-xr+(xr-xl)/100
		dx<-(xmax-xmin)/(knots.no-1)
		knots<-seq(xmin-spline.degree*dx,xmax+spline.degree*dx,by=dx)
		call$knots <- knots
	}
	
	
	B<-spline.des(knots,x,spline.degree+1)$design
	
	
	#generate Penalization-Matrix
	D<-diag(k)
	if((d<-min(diff.ord,spline.degree))<diff.ord) warning(paste("order of differences > degree of splines:\n new order of differences=",d,"\n"))
	if(d > 0) for(i in 1:d) D<-diff(D)
	call$diff.ord <- d
	
	#reparametrization: X = unpenalized part, Z =penalized part
	X <-rep(1,k)
	if(diff.ord>1) {for(i in 2:diff.ord) X <-cbind(X, knots[1:k]^(i-1)) }
	
	X <- (B%*%X)[,-1, drop=FALSE]
	Z  <- B%*%t(D)%*%solve(tcrossprod(D))
	
	res <- list(X=unname(X), Z=unname(Z))
	attr(res, "call") <- as.call(call)	
	return(res)
}
