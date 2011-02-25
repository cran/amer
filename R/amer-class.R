setClass("amer", representation = representation(smooths = "list"),
		contains = "mer")

setMethod("summary", signature(object = "amer"),
		function(object, ...)
		{
			REML <- object@dims["REML"]
			fcoef <- fixef(object)
			vcov <- vcov(object)
			corF <- vcov@factors$correlation
			dims <- object@dims
			coefs <- cbind("Estimate" = fcoef, "Std. Error" = corF@sd) #, DF = DF)
			llik <- logLik(object, REML)
			dev <- object@deviance
			mType <- if((non <- as.logical(length(object@V)))) "NMM" else "LMM"
			if (gen <- as.logical(length(object@muEta)))
				mType <- paste("G", mType, sep = '')
			mName <- switch(mType, LMM = "Additive", NMM = "Nonlinear additive",
					GLMM = "Generalized additive",
					GNMM = "Generalized nonlinear additive")
			method <- {
				if (mType == "LMM")
					if(REML) "REML" else "maximum likelihood"
				else
					paste("the", if(dims["nAGQ"] == 1) "Laplace" else
										"adaptive Gaussian Hermite",
							"approximation")
			}
			AICframe <- data.frame(AIC = AIC(llik), BIC = BIC(llik),
					logLik = as.vector(llik),
					deviance = dev["ML"],
					REMLdev = dev["REML"],
					row.names = "")
			if (is.na(AICframe$REMLdev)) AICframe$REMLdev <- NULL
			varcor <- lme4:::VarCorr(object)
			REmat <- lme4:::formatVC(varcor)
			if (is.na(attr(varcor, "sc")))
				REmat <- REmat[-nrow(REmat), , drop = FALSE]
			
			if (nrow(coefs) > 0) {
				if (!dims["useSc"]) {
					coefs <- coefs[, 1:2, drop = FALSE]
					stat <- coefs[,1]/coefs[,2]
					pval <- 2*pnorm(abs(stat), lower = FALSE)
					coefs <- cbind(coefs, "z value" = stat, "Pr(>|z|)" = pval)
				} else {
					stat <- coefs[,1]/coefs[,2]
					##pval <- 2*pt(abs(stat), coefs[,3], lower = FALSE)
					coefs <- cbind(coefs, "t value" = stat) #, "Pr(>|t|)" = pval)
				}
			} ## else : append columns to 0-row matrix ...
			new("summary.mer",
					object,
					methTitle = paste(mName, "mixed model fit by", method),
					logLik = llik,
					ngrps = sapply(object@flist, function(x) length(levels(x))),
					sigma = lme4:::sigma(object),
					coefs = coefs,
					vcov = vcov,
					REmat = REmat,
					AICtab= AICframe
			)
		})## summary()


#if(FALSE){
#	
	setGeneric("predict",
			function(object, ...)
				standardGeneric("predict")
	)
	## returns newdata with additional columns for:
	##  - the components of the linear predictor X*beta, named as in fixef(object)
	##  - the components of the random effects Z_i%*%b_i, named after the group levels
	##  - the smooths, named as in names(objects@smooths)  
	##  - column "fit" contains:
	##		if type="response": the fitted values on the scale of the response;
	##		if type="linpred": sum of X%*%beta + Z %*% b + \sum f_s
	## if newdata contains additional levels of grouping factors not present in the original model, 
	##	these are assigned random effects of zero.
	setMethod("predict", signature(object = "amer"),
			function(object, newdata, type=c("response", "linpred", "terms"), ...)
			{
				
				type <- match.arg(type)
				objcall <- object@call
				
				
				#add column for reponse if missing and pad factors s.t. they have more than one level
				# (lmerSetup aborts if not)
				padding <- NULL
				addResp <- FALSE 
				if(!((response <- deparse(object@call$formula[[2]])) %in% colnames(newdata))){
					args <-list(newdata, rep(c(0,1), l=nrow(newdata)))
					names(args) <- c("", response)
					newdata <- do.call(cbind, args)
					addResp <- TRUE
				} 
				if(any(isFactor <- sapply(as.list(newdata),is.factor))){
					if(any(constant <- 
									sapply(as.list(newdata[, isFactor, drop=FALSE]), 
											function(x){
												length(unique(x)) == 1
											}))){
						padding <- newdata[1,]
						for(f in which(isFactor)[constant]){
							lvls <- levels(newdata[, f])
							use <- lvls[lvls != unique(newdata[, f])][1]
							padding[,f] <- use
						}
						newdata <- rbind(newdata, padding)
					}
				}
				objcall$data <-  newdata
				
				newObj <- do.call("amerSetup", as.list(objcall)[-1])$m

				#add columns for dropped levels & reorder
				if(any(dropped <- !(names(object@fixef) %in% names(newObj$fr$fixef)))){
					padCols <- matrix(0, nrow=nrow(newdata), ncol = sum(dropped))
					newObj$fr$X <- cbind(newObj$fr$X, padCols)
				}
				colnames(newObj$fr$X) <- c(names(newObj$fr$fixef), names(object@fixef)[dropped])
				newObj$fr$X <- newObj$fr$X[, names(fixef(object)), drop=FALSE]
				
				#find indices for fixed/random/smooth effects
				indFixed <- (1:ncol(newObj$fr$X))[-unlist(sapply(object@smooths, function(sm) {
											return(
													ifelse(eval(sm$allPen),
															ncol(newObj$fr$X)+1,
															attr(sm, "indUnpen"))
											)
										}))] 
				indSmoo <- unlist(sapply(object@smooths, function(sm) attr(sm, "indGrp")))  
				indRan <- (1:length(newObj$FL$fl))[!(names(newObj$FL$fl) %in% names(object@flist)[indSmoo])] 
				
				#X*beta
				if(length(indFixed)){
					fixed <- t(t(newObj$fr$X[,indFixed, drop=F]) * fixef(object)[indFixed, drop=F])
					colnames(fixed) <- paste("lp.",names(fixef(object)[indFixed]),sep="")
				} else fixed <- NULL 
				
				
				#Z%*%b
				if(length(indRan)){
					random <- sapply(indRan, function(i){
								thisName <- names(newObj$FL$fl)[attr(newObj$FL$fl, "assign")[i]]
								thisInd <- which(names(newObj$FL$fl) == thisName)
								useZRows <- (newObj$FL$fl[[thisName]] %in% levels(object@flist[[thisName]]))
								useZCols <- (levels(newObj$FL$fl[[thisName]]) %in% levels(object@flist[[thisName]]))
								useBCols <- (levels(object@flist[[thisName]]) %in% levels(newObj$FL$fl[[thisName]]))
								##FIXME: this works only if there are either: 
								## - levels of grouping in newdata that are not in orig data
								## - levels of grouping in orig data that are not in newdata
								## - neither,  but NOT for both...
								tmp <- rep(0, nrow(newdata))
								tmp[useZRows] <- as.numeric(t(newObj$FL$trms[[thisInd]]$Zt)[useZRows, useZCols, drop=FALSE] %*% 
												unlist(ranef(object)[[thisName]])[useBCols, drop=FALSE])
								tmp
							}) 
					colnames(random) <- paste("blup.",names(newObj$FL$fl)[indRan], sep="")
				} else random <- NULL
				

				fctList <- getF(object, newdata=newdata, addConst=FALSE)
				
				smooth <- sapply(1:length(fctList), 
						function(i){
							if(length(fctList[[i]]) == 1) return(fctList[[i]][[1]]$fhat)
							else { #separate by level of by-variable
								ff <- sapply(fctList[[i]], function(x) x$fhat)
								grp <- fctList[[i]][[1]][, deparse(object@smooths[[i]]$by)]
								levels(grp) <- 1:nlevels(grp)
								return(ff[cbind(1:length(grp),grp)])
							}
						} ) 			
				colnames(smooth) <- names(fctList) 
				
				fits <- do.call(cbind, list(fixed, random, smooth))
				
				res <- cbind(newdata, fits)
				if(addResp) res <- res[, colnames(res) != response]
				if(!is.null(padding)) {
					res <- res[-nrow(res), ]
					fits <- fits[-nrow(fits), ]
				}	
				
				if(type=="terms") return(res)
				else {
					pred <- rowSums(fits)
					if(type=="response" && !is.null(object@call$family)) {
						linkinv <- eval(parse(text = paste(object@call$family,"()", sep = "")))$linkinv
						pred <- do.call(linkinv, list(eta = pred))
					}
					return(cbind(res, fit = pred))
				}
			}
	)## predict() 
#}
