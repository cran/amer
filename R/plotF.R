#' plot estimated function values + pointwise confidence intervals for an amer-Fit 	
#'
#' @param object a fitted additive (mixed) model of class \code{\link{amer-class}} 
#' @param which (optional) an integer vector or a character vector of names giving the smooths for which plots are desired. Defaults to all. 
#' @param n fitted values for a regular grid with n values in the range of the respective covariates are calculated   
#' @param interval what mehod should be used to compute pointwise confidence/HPD intervals. See \code{\link{getF}}.  
#' @param addConst boolean should the global intercept and intercepts for the levels of the by-variable be included in the fitted values (and their CIs) 
#' @param trans a function that should be applied to the fitted values and ci's before plotting (e.g. the inverse link function to get plots on the scale of the reponse)
#' @param level level for the confidence/HPD intervals
#' @param sims how many iterates should  be generated for the MCMC-based HPD-intervals
#' @param auto.layout automagically set plot layout via \code{par()$mfrow}
#' @param rug add \code{\link{rug}}-plots of the observed covariate locations
#' @param legendPos a (vector of) keyword(s) where to put labels of by-variables (see \code{\link[graphics]{legend}}). "none" if you don't want a legend.  
#' @param ... arguments passed on to the low-level plot functions (\code{plot}, \code{matlines}), \code{legend}, and \code{title}
#' @return invisibly returns the result of the \code{\link{getF}}-call used to do the plots
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{getF}}; \code{\link{plotF}},  \code{tests/optionsTests.r} and the vignette for examples 

plotF <- function(object, which, n=100, interval = "RW", addConst = TRUE, trans=I,  
		level = 0.9, sims = 1000, auto.layout = TRUE, rug = TRUE, legendPos="topright", ...)
{
# FIXME: add option for centering function estimates at zero/ anchoring at mean of all other covariates?
# FIXME: valid confints for trans? 	
	terms <- object@smooths
	if(missing(which)) which <- seq_along(terms)
	if(is.character(which)) {
		which <- match(which, names(terms))
		if(any(nas <- is.na(which))) 
			warning("entry ", paste(which[nas], collapse=", "), " in 'which' did not match any function names in ", safeDeparse(object) ,".")
		which <- which[!nas]
	}
	if(length(legendPos) != length(which)) legendPos <- rep(legendPos, length=length(which))
	if(length(addConst) != length(which)) addConst <- rep(addConst, length=length(which))
	
	interval <- toupper(interval)
	
	res <- getF(object = object, which = which, n = n, newdata=NULL, interval = interval, addConst = addConst, level = level, sims = sims)
	allPen <- sapply(object@smooths[which], function(x) eval(x$allPen))
	
	if(auto.layout){
		oldpar <- NULL
		on.exit(par(oldpar))
		nf <- length(res)
		oldpar <- par(mfrow = set.mfrow(Nparms=nf))
	}
	
	
	plot1F <- function(x, interval, legendPos, ...){
		dots <- list(...)
		fhat <- sapply(x, function(x){trans(x$fhat)})
		if(interval) ci <- lapply(x, function(x){trans(cbind(x$lo, x$hi))})
		cov <- sapply(x, function(x){x[,1]})
		if(is.null(dots$ylim)){
			ylim <- range(fhat)
			if(interval) ylim <- range(ylim, ci)
		} else {
			ylim <- dots$ylim
			dots <- dots[names(dots)!="ylim"]
		}	
		if(is.null(dots$xlim)){
			xlim <- range(cov)
		} else {
			xlim <- dots$xlim
			dots <- dots[names(dots)!="xlim"]
		}
		do.call(plot, c(list(x=-2*xlim[1]-10, y=0, ylim = ylim, xlim = xlim, ylab="", xlab=colnames(x[[1]])[1]), dots))
		do.call(matlines, c(list(x = cov, y = fhat, col=1:length(x), lty=1), dots))
		if(interval) for(i in 1:length(ci)) do.call(matlines, c(list(x=cov[,i], y=ci[[i]], col=i, lty=2),dots))
		if(legendPos != "none" && length(x)>1){
			do.call(legend,c(list(x=legendPos, legend=names(x), lty=1, col=1:length(x)),dots))
		}	
	}
	dots <- list(...)
	for(i in seq_along(res)){
		#FIXME: change as CIs for allPen become available
		plot1F(res[[i]], interval = !(interval=="NONE")&&!allPen[i], legendPos = legendPos[i], ...)
		if(is.null(dots$ylab)){
			ylab <- ifelse(addConst[i], paste(names(res)[i], "+ const"), names(res)[i])
			if(any(trans(-2:2)!= (-2:2))) ylab <- paste(safeDeparse(match.call()$trans),"(",ylab,")",sep="")
		} else{
			ylab <-dots$ylab
			dots <- dots[names(dots)!="ylab"]
		}	
		do.call(title, c(list(ylab= ylab), dots))
		if(rug){
			if(length(res[[i]])==1){
				rug(object@frame[,colnames(res[[i]][[1]])[1]], ...)
			} else {
				nlvls <- length(res[[i]]) 
				lvls <- levels(object@frame[, safeDeparse(object@smooths[[i]]$by)])
				for(j in 1:nlvls){
					use <- object@frame[,safeDeparse(object@smooths[[i]]$by)] == lvls[j]
					rug(object@frame[use, colnames(res[[i]][[1]])[1]], col =j, ...)
				}	
			}	
		} 
	}
	invisible(res)
}
# plotF(test1, trans=exp)
