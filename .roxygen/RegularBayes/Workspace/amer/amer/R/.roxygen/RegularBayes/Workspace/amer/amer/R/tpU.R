#' Generate a modified truncated power basis for penalized spline smoothing.
#'
#' @param x covariate for the smooth function
#' @param degree integer: degree of truncated polynomials (0: piecewise constant, 1: piecewise linear etc..)
#' @param k integer: dimensionality of the basis (i.e.: number of knots + degree)
#' @param unpen integer: degree of the unpenalized nullspace, must be lower than degree: 1 for pen. deviations from linearity, 2 for pen. deviations from quadratic from etc.
#' @param by factor variable: estimate separate functions for each level - this assumes standard treatment contrasts for the supplied factor.
#' @param allPen boolean: if TRUE, make design for group-specific curves with common smoothing parameter: all parameters (including the normally unpenalized basis functions in X) are penalized, every level of "by" has the same amount of smoothing. If FALSE, make design for separate curves for each by-level: separate smoothing parameters for every level of "by", unpenalized estimates for the coefficients associated with X
#' @param varying numeric: if not NULL, a varying coefficient model is fit: f(x,varying) = f(x)*varying
#' @param diag logical: force a diagonal covariance-matrix for the random effects for X if \code{allPen=TRUE}? 
#' @param knots vector of knot locations (optional). Defaults to quantile-based knot placement at the \eqn{(i+1)/(k+2)}-quantiles for \eqn{i=1,\dots,k}. 
#' @param centerscale numeric(2): center&scale x by these values if not NULL
#' @param scaledknots boolean:	are knots given for the rescaled x-values?
#' @note  This is a more detailed implementation of the example on how to define additional basis generating functions in the vignette.
#' @return list with entries:
#'				\code{"X"}: \code{n x unpen} design matrix for unpenalized part (without intercept)  
#' 				\code{"Z"}: \code{n x (k-unpen)} design matrix for penalized part 
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{tp}}
tpU <-
		function(x, degree=2, k = 15, unpen=1, by=NULL, allPen = FALSE, varying = NULL, diag=FALSE,
				knots= quantile(x, probs = (2:(k - degree + 1))/(k - degree  + 2)), centerscale=NULL, scaledknots=FALSE)
{
	call <- as.list(expand.call())
	
	stopifnot(is.numeric(x), is.factor(by)||is.null(by), is.numeric(varying)||is.null(varying), degree >= 0, unpen <= degree)
	
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
	if((knots[1]<min(x)||(knots[k-degree]>max(x)))) warning("some knots outside range of variable.")
	
	if(is.null(by) && allPen) stop("allPen = TRUE only makes sense for smooths with a by-variable.")
	
	#design for unpenalised part: global polynomial trends (no intercept)
	if(degree>0){	
		X <- outer(x, 1:degree, "^")#poly(x, degree)
	} else{
		X <- matrix(nrow=length(x), ncol=0)
	}	
	#design for penalised part: 
	Z <- outer(x,knots,"-")^degree*outer(x,knots,">")
	
	# adapt design for unpen option
	if(unpen != degree){
		Xu <- X[,(1:degree) > unpen, drop=F]
		X <- X[,(1:degree) <= unpen, drop=F]
		Z <- cbind(Xu, Z)
	}
	
	
	res <- list(X=unname(X), Z=unname(Z), knots=knots)
	attr(res, "call") <- as.call(call)
	
	return(res)
}


