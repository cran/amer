#' Generate a truncated power basis for penalized spline smoothing.
#' 
#' Truncated power bases have \code{degree} unpenalized basis functions, namely \eqn{x^1,\dots, x^{degree}} and \eqn{k-}\code{degree} penalized basis functions that contain the positive part \eqn{(x-\kappa_j)^{degree}} for knots \eqn{\kappa_j, j=1,dots,k-}\code{degree}.
#' This function can be used as a reference when implementing other \code{basisGenerators} that can be used for \code{\link{amer}}-fits. 
#' All such functions need to return a list of at least X (unpenalized basis functions, a matrix with zero columns if there are none of those), and Z (penalized basis functions)
#' that has a \code{call}-attribute with the expanded call returned by \code{\link{expand.call}()}. All such functions need to have at least arguments \code{x, by, allPen, diag} and \code{varying}.
#' See also section 4.4 in the vignette for an example on how to write your own basis-generating functions.
#'
#' @param x covariate for the smooth function
#' @param degree integer: degree of truncated polynomials (0: piecewise constant, 1: piecewise linear etc..)
#' @param k integer: dimensionality of the basis (i.e.: number of knots + degree)
#' @param by factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.
#' @param allPen boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing 
#'                         if FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X
#' @param varying numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying
#' @param diag logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}? 
#' @param knots vector of knot locations (optional). Defaults quantile-based knots at the \eqn{(i+1)/(k+2)}-quantiles 
#' 		  for \eqn{i=1,\dots,k}.
#' @param centerscale numeric(2): center&scale x by these values if not NULL
#' @param scaledknots boolean:	are knot locations given for the rescaled x-values?
#' @return list with entries:
#'			\code{"X"}: \code{n x degree} design matrix for unpenalized part (without intercept) (or a list of those for every level of by if allPen=F), 
#' 			\code{"Z"}: \code{n x (k-degree)} design matrix for penalized part (or a list of those for every level of by if allPen=F),
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{tp}}
tp <-
function(x, degree=1, k = 15, by=NULL, allPen = FALSE, varying = NULL, diag=FALSE,
		knots= quantile(x, probs = (2:(k - degree + 1))/(k - degree  + 2)), centerscale=NULL, scaledknots=FALSE)
{
	call <- as.list(expand.call())
	
	stopifnot(is.numeric(x), is.factor(by)||is.null(by), is.numeric(varying)||is.null(varying), degree >= 0)
	
	degree <- as.integer(degree); call$degree <- degree
	
	knots <- eval(knots)
	if(is.null(centerscale)){
		x <- scale(x)
		#make sure prediction data uses the same center/scale as the data used to fit the model:
		call$centerscale <- c(attr(x, "scaled:center"),attr(x, "scaled:scale"))
		x <- as.vector(x)
	} else x <- (x - centerscale[1])/centerscale[2]
	if(!scaledknots){
		knots <- (knots - call$centerscale[1])/call$centerscale[2]
		call$scaledknots <- TRUE
	}	
	
	if(length(unique(knots))!=length(knots)) warning("duplicate knots detected and removed.")
	knots <- sort(unique(knots))

	call$knots <- knots
	if(k != length(knots)+ degree){
		k <- length(knots) + degree; call$k <- k
		warning("set k to ", k," to conform with given knots and degree.")
	}
	if((knots[1]<min(x)||(knots[k-degree]>max(x)))) warning("knots outside range of variable.")
	
	if(is.null(by) && allPen) stop("allPen = TRUE only makes sense for smooths with a by-variable.")
	
	#design for unpenalised part: global polynomial trends (no intercept)
	if(degree>0){	
		X <- outer(x, 1:degree, "^")#poly(x, degree)
		#colnames(X) <- paste(as.character(call$x),".fx",1:NCOL(X),sep="")
	} else{
		X <- matrix(nrow=length(x), ncol=0)
	}	
	#design for penalised part: 
	Z <- outer(x,knots,"-")^degree*outer(x,knots,">")
	
	
	res <- list(X=X, Z=Z, knots=knots)
	attr(res, "call") <- as.call(call)

	return(res)
}

