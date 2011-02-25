amerSetup <-
function(formula, data, family, REML, control = list(), 
		start, verbose, 
#		subset, weights, 
#		na.action, offset, 
		contrasts, 
		basisGenerators, bySetToZero = T){
#set up unfitted amer-object, similar to lme4:::lmerFrames:
#1. create spline bases [tp(), expandBasis()]
#2. manipulate formula for lmer [subFcts()]
#3. augment model frame with unpenalized basis functions and 
#	fake grouping factors as placeholders for penalized basis fcts [expandMf()] 
#4. set up model with augmented model frame 
#5. overwrite design matrices for fake random intercepts with designs for the penalized basis functions	[subAZ()]
# returns list: 
# 	m: unfitted model object, 
#	fct: list of spline bases (see tp(), expandBasis()) 	
#	fctterm: list of smooth term specifications

	call <- match.call()
	
	formula <- eval(call$formula)
	
	
	# identify smooth terms
	tf <- terms.formula(formula, specials = eval(call$basisGenerators, parent.frame(2)))
	f.ind <- unlist(attr(tf,"specials"))
	n.f <- length(f.ind)
	
	if(n.f){
		fixed.formula <- as.character(lme4:::nobars(formula)[3])
		
		if(!is.null(bars <- lme4:::findbars(formula))){
			random.formula <- paste("(", bars, ")", collapse=" + ")
			rhs <- paste(fixed.formula, random.formula, sep=" + ")
		} else rhs <- fixed.formula
		
		fctterm <- fct <- vector(mode="list", length= n.f)
		
		# extract function terms from formula ...
		for(i in 1:n.f) fctterm[[i]] <- attr(tf,"variables")[[f.ind[i]+1]]
		
		# make sure data is available
		if(is.null(call$data)){
			myGsub <- function(x, pattern, replacement) gsub(as.character(pattern), replacement, x)
			fake.formula <- paste(formula[[2]], "~", paste(lapply(as.character(bars), myGsub, pattern="\\|", replacement="+"), collapse=" + ")) 
			getSmoothCovariates <- function(term) as.call(term)[[2]]
			fake.formula <- formula(paste(fake.formula, "+", paste(lapply(fctterm, getSmoothCovariates), collapse = " + ")))
			data <- model.frame(fake.formula)
		} else data <- eval.parent(call$data, n=2)
		
		
		# ... set up designs (see: tp())...
		#browser()
		fct <- lapply(fctterm, eval, envir = data, enclos=parent.frame(2))
		for(i in seq_along(fct)) fct[[i]] <- expandBasis(fct[[i]], 
					eval(attr(fct[[i]], "call")$by, data),
					eval(attr(fct[[i]], "call")$varying, data),
					bySetToZero) 
			
		
		#naming scheme: "f.x" if allPen=F, else "f.x.by" OR "f.xXvarying", "f.xXvarying.by" respectively
		names(fct) <- names(fctterm) <- 
				paste("f.", 
						lapply(fct, function(x){
									paste(as.character(attr(x,"call")$x),
											ifelse(!is.null(eval(attr(x,"call")$varying, data)),  paste("X", deparse(attr(x, "call")$varying), sep=""), ""),
											ifelse(eval(attr(x,"call")$allPen),  paste(".", deparse(attr(x, "call")$by), sep=""), ""),
											sep="")}),
						sep="")						
		
		#... adjust formula ...
		rhs <- subFcts(rhs, fctterm, fct, data)
		
		#... add unpenalized part of bases and fakeFactors to model frame
		data <- expandMf(data, fct)
	} else stop("no smooth terms given - use lmer instead.")	
	
	#set up lmer-call with fake formula and augmented data
	call[[1]] <- as.name("lmer")
	call$doFit <- FALSE
	call$data <- as.name("data")
	call$formula <- as.formula(paste(formula[[2]],"~",rhs))
	#FIXME: should intercept be enforced? - may not always make sense, esp. if by-variable is given (also need to change in lamer if changed here)
	#### make sure there is an intercept: 
	###call$formula <- update(call$formula, .~.+1) 
	m <- eval(call, data)
	#make sure none of the variables are dropped
	m$fr$mf <- data 
	
	#replace A, Z, Zt of fake factors with the penalized bases for the smooth terms
	m <- subAZ(m, fct)
	# make fctterm more informative: use expanded calls from expand.call()
	fctterm <- lapply(fct, function(x) attr(x, "call"))
	
	return(list(m=m, fct=fct, fctterm=fctterm))
}

