#' get the estimated function values from an amer-Fit 	
#'
#' @param object a fitted additive (mixed) model of class \code{\link{amer-class}} 
#' @param which (optional) an integer vector or a character vector of names giving the smooths for which fitted values are desired. Defaults to all. 
#' @param n if no \code{newdata} is given, fitted values for a regular grid with n values in the range of the respective covariates are returned   
#' @param newdata An optional data frame in which to look for variables with which to predict
#' @param interval what mehod should be used to compute pointwise confidence/HPD intervals: RW= bias-adjusted empirical bayes, MCMC uses \code{\link[lme4]{mcmcsamp}} 
#' @param addConst boolean should the global intercept and intercepts for the levels of the by-variable be included in the fitted values (and their CIs) can also be a vector of the same length as \code{which} 
#' @param varying value of the\code{varying}-covariate (see \code{\link{tp}}) to be used if no newdata is supplied. 
#' 			Defaults to 1. 
#' @param level level for the confidence/HPD intervals
#' @param sims how many iterates should  be generated for the MCMC-based HPD-intervals
#' @note The formula used for the pointwise bias-adjusted CIs is taken from Ruppert and Wand's  'Semiparametric Regression' (2003), p. 140. 
#' 			These leave out the uncertainty associated with the variance component estimates. 
#' 			MCMC-intervals based on results from \code{\link[lme4]{mcmcsamp}} don't seem to be very reliable yet and should be used with caution, especially for more complex models.
#' 			  
#' @return a list with one \code{data.frame} for each function, giving \code{newdata} or the values of the generated grid plus the fitted values (and confidence/HPD intervals)
#' 		if MCMC-intervals were rquested, the list has an attribute "mcmc" containing the result of the call to \code{\link[lme4]{mcmcsamp}}, a \code{\link[lme4]{merMCMC-class}} object.	
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link{plotF}},  \code{tests/optionsTests.r} and the vignette for examples 
getF <- function (object, which, n=100, newdata=NULL, interval = c("NONE", "MCMC", "RW"), addConst = TRUE, varying=1, 
		level = 0.9, sims = 1000)
{

	stopifnot(class(object)=="amer", is.null(newdata)||is.data.frame(newdata), n>0, sims>0, level < 1, level > 0)
	
	
	terms <- object@smooths
	interval <- toupper(interval)
	interval <- match.arg(interval)
	
	n <- as.integer(n)
	sims <- as.integer(sims)
	if(missing(which)) which <- seq_along(terms)
	if(is.character(which)) {
		which <- match(which, names(terms))
		if(any(nas <- is.na(which))) 
			warning("entry ", paste(which[nas], collapse=", "), " in 'which' did not match any function names in ", deparse(object) ,".")
		which <- which[!nas]
	}
	if(length(addConst) != length(which)) addConst <- rep(addConst, length=length(which))
	
	if(interval=="MCMC"){
		#FIXME: these are sometimes HUGE compared to RW, with strange shapes and/or not overlapping the posterior means/modes 
		#		- cannot be just because variability of variance estimates enters here but not in RW (?)
		#		- does mcmcsamp break for data that's not really grouped?
		warning("mcmcsamp() may not work reliably yet -- check the traceplots!.")
		cat("starting", sims,"MCMC iterations for posterior intervals:\n")
		mcmc <- mcmcsamp(object, n=sims, saveb=T)
		cat("		... done.\n")
		
		ci.MCMC<- function(base, terms, i, j, mcmc, level, addConst){

			indUnpen <- if(addConst){
						c(attr(terms[[i]],"indConst")[[j]], attr(terms[[i]],"indUnpen")[[j]])
					} else attr(terms[[i]],"indUnpen")[[j]]	
			fhat <- base$X[[j]] %*% mcmc@fixef[indUnpen,] +
					base$Z[[j]] %*% mcmc@ranef[attr(terms[[i]],"indPen")[[j]],]
			ci <- HPDinterval(fhat, prob=level)
			#probs <- c((1-level)/2, level + (1-level)/2)
			#ci <- t(apply(fhat, 1, quantile, probs = probs, na.rm=T)
			colnames(ci) <- c("lo", "hi")
			return(ci)
		}
	}
	
	ci.RW <- function(fhat, base, terms, i, j, object, level, addConst){
		fctV <-
				function(m, indGrp, indPen, indUnpen){
			#Cov(hat.beta, hat.b-b) for bias-adjusted empirical Bayes CIs (s. Ruppert/Wand(2003), Semiparametric Regression, p. 138 f.):
			#use V = cov(hat.fixef, hat.ranef) = sigma.eps^2 (C'C + sigma.eps^2/sigma.b^2 D)^-1; C=[XZ]*sqrt(W), D = blockdiag(0, I_dim(b))
			
			C <- cBind(m@X[,indUnpen, drop=F], t(m@Zt[indPen,]))
			if(length(m@var)) C <- C * sqrt(1/m@var)
			
			V <- crossprod(C)
			# FIXME: this works only for scalar ST and length(indGrp=1)- don't use if allPen=T! 
			if(m@ST[[indGrp]] > 0){
				diag(V)[-(1:length(indUnpen))] <- diag(V)[-(1:length(indUnpen))] + m@ST[[indGrp]]^-2
			} else {
				#FIXME: HACK: V is not invertible for var(ranef)=0, use var(ranef)=10^-9 instead 
				diag(V)[-(1:length(indUnpen))] <- diag(V)[-(1:length(indUnpen))] + 10^9
			}
			return(lme4:::sigma(m)^2 * solve(V))
		}
		
		z <- qnorm(1-(1-level)/2)
		indUnpen <- if(addConst){
					c(attr(terms[[i]],"indConst")[[j]], attr(terms[[i]],"indUnpen")[[j]])
					} else attr(terms[[i]],"indUnpen")[[j]]	
		cV <- chol(fctV(object, attr(terms[[i]],"indGrp")[[j]], 
						attr(terms[[i]],"indPen")[[j]], indUnpen))
		
		C <- as(cBind(base$X[[j]], base$Z[[j]]), "sparseMatrix")
		sd <- apply(C, 1, function(x, cV){
					return(sqrt((crossprod(cV%*%x))@x))	
				}, cV=cV)
		
		ci <- cbind(fhat - z*sd, fhat + z*sd)
		colnames(ci) <- c("lo", "hi")
		return(ci)
	}
	
	ans <- vector(mode="list", length=length(which))
	if(interval=="MCMC") attr(ans, "mcmc") <- mcmc 
	names(ans) <- names(terms)[which]
	indWhich <- 1
	
	for(i in which){
		#################################
		# set up / check newdata
		################################
		if(!is.null(terms[[i]]$by)){
			lvls <- levels(object@frame[, deparse(terms[[i]]$by)])#FIXME: in amerSetup: this will fail if terms[[i]]$x was only in the workspace but not in the supplied data.frame for the original call.
			hasBy <- TRUE
		} else hasBy <-FALSE	
		hasVarying <- !is.null(terms[[i]]$varying)
		
		if(is.null(newdata)){
			grid <- TRUE
			#FIXME: adapt this for 2d/3d-smooths!!
			#get range of covariates + sequence of values
			lim <- range(object@frame[, deparse(terms[[i]]$x)], na.rm=T) #FIXME: in amerSetup: this will fail if terms[[i]]$x was only in the workspace but not in the supplied data.frame for the original call.
			newX <- seq(lim[1], lim[2], l=n)
			if(hasBy){
				newBy <- factor(rep(lvls, length=n), labels=lvls)
				data <- data.frame(newX, newBy)
				colnames(data) <- c(deparse(terms[[i]]$x), deparse(terms[[i]]$by))
			} else {
				data <- data.frame(newX)
				colnames(data) <- deparse(terms[[i]]$x)
			}
			if(hasVarying){
				#varying covariate is set value of varying
				data <- cbind(data, rep(varying, nrow(data)))
				colnames(data)[NCOL(data)] <- deparse(terms[[i]]$varying)
			}
		} else {
			grid <- FALSE
			data <- newdata
			n <- nrow(newdata)
			vnames <- deparse(terms[[i]]$x)
			if(hasBy) vnames <- c(vnames,deparse(terms[[i]]$by))
			if(hasVarying) vnames <- c(vnames, deparse(terms[[i]]$varying))
			if(any(nas <- is.na(match(vnames, colnames(data))))) 
				stop("variable ", paste(vnames[nas], collapse=", "), "not found in given data.")
			data <- data[, colnames(data) %in% vnames, drop=F]
		}
		
		#################################
		#create basis:
		#################################
		base <- eval(terms[[i]], data)
		base <- expandBasis(base, 
				by = eval(attr(base, "call")$by, data),  
				varying = eval(attr(base, "call")$varying, data),
				bySetToZero = !grid)

		# need to modify Z, X for allPen-Fits
		if(eval(terms[[i]]$allPen)){
			nlvl <- length(lvls)
			base0 <- base
			#where are the random effects for the penalized spline functions:
			indZ <- lme4:::reinds(object@Gp)[[attr(terms[[i]],"indGrp")[[1]][1]]] ##can use the first because the spline will have more levels (grps*(p-d)) than the grouping factor (grps) 
			#how many penalized spline functions per level of by
			dimOneZ <- length(indZ)/length(lvls)
			useZ <- 1:dimOneZ
		}	
		
		nf <- ifelse(hasBy, length(lvls), 1)
		ans[[indWhich]] <- vector(mode="list", length = nf)
		names(ans[[indWhich]]) <-  if(hasBy) paste(deparse(terms[[i]]$by), lvls, sep="") else names(base$X)
		for(j in seq_along(ans[[indWhich]])){
		#################################
		#calculate fits and cis
		#################################	

			if(eval(terms[[i]]$allPen)){
				ansInd <- 1
				#need to:
				#-append to Z extra columns for the random effects for base$X	 (+ a random intercept for the by-levels)
				#-set columns in base$Z not relevant for the current by-level to 0
				#-set base$X to zero
				base$Z[[1]] <- base0$Z[[1]]
				base$Z[[1]][,-useZ] <- 0
				#step to next block:
				useZ <- useZ + dimOneZ
				lvlInd <- rep(0, nlvl)
				lvlInd[j] <- 1
				
				X <- cBind("(Intercept)"=1, base0$X[[1]])
				if(!grid){
					use <- eval(terms[[i]]$by, data)==lvls[j]
					X[!use,] <- 0
				} 
				if(eval(terms[[i]]$diag)){
					#lmer switches order of terms in X, need to permuite X-columns accordingly: 
					uNames <-  unlist(unique(sapply(object@ST[attr(terms[[i]],"indGrp")[[1]][-1]], dimnames)))
					X <- X[,uNames]
				}
				
				base$Z[[1]] <-  as(cBind(base$Z[[1]], kronecker(X, t(lvlInd), FUN = "*")),"sparseMatrix")
				base$X[[1]] <- matrix(0, nrow=n, ncol=0)
			} else {
				ansInd <- j
				if(!grid && hasBy){
						#remove unnecessary rows from design
						use <- eval(terms[[i]]$by, data)==lvls[j]
						base$X[[ansInd]] <- base$X[[ansInd]][use,]
						base$Z[[ansInd]] <- base$Z[[ansInd]][use,]
				}
			}	
		
			if(addConst[i]){
			#add columns for constant terms to X, append indUnpen:	
				byColumn <-if(hasBy && paste(deparse(terms[[i]]$by),lvls[j],sep="") %in% names(object@fixef)){
				 				 rep(1, nrow(base$X[[ansInd]]))
							} else numeric(0) 
				base$X[[ansInd]] <- cBind(byColumn, base$X[[ansInd]])	
				if("(Intercept)" %in% names(object@fixef)[attr(terms[[i]],"indConst")[[ansInd]]])
					base$X[[ansInd]] <- cBind(1,base$X[[ansInd]])
				
				indUnpen <- c(attr(terms[[i]],"indConst")[[ansInd]], attr(terms[[i]],"indUnpen")[[ansInd]])
			} else indUnpen <- attr(terms[[i]],"indUnpen")[[ansInd]]	
		
			fhat <- base$X[[ansInd]] %*% 
					object@fixef[indUnpen, drop = F] + 
					base$Z[[ansInd]] %*% 
					object@ranef[attr(terms[[i]],"indPen")[[ansInd]], drop = F] 
			
			if(!eval(terms[[i]]$allPen)){ 
				ci <- switch(interval,
						"NONE" = matrix(NA, nrow=nrow(base$X[[ansInd]]), ncol=0),
						"RW" = ci.RW(fhat, base, terms, i, j, object, level, addConst[i]),
						"MCMC" = ci.MCMC(base, terms, i, j, mcmc, level, addConst[i]))
			} else {
				#TODO: implement CIs for fits with allPen = T
				if(interval!="NONE" && j==1) warning("CIs for fits with allPen = T not yet implemented.")
				ci <- matrix(NA, nrow=nrow(base$X[[ansInd]]), ncol=0)
			}
			dataJ <- if(grid){
				data[,!(colnames(data)==deparse(terms[[i]]$by)), drop=F]
			} else {
				if(!grid && hasBy){
					data[use,] 
				} else data  
			}	
			ans[[indWhich]][[j]] <- data.frame(dataJ, fhat= as.matrix(fhat), ci)
		}# end for j
		indWhich <- indWhich + 1
	}#end for i
	return(ans)
}







