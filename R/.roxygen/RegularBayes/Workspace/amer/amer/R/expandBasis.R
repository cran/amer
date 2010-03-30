expandBasis <-
function(basis, by, varying, bySetToZero = T){
#multiply X, Z with varying
#set up lists of design matrices if by-variable and/or allPen are given: 	
	allPen <- eval(attr(basis, "call")$allPen, parent.frame())
	X <- basis$X
	Z <- basis$Z
	
	
	xName <- deparse(attr(basis, "call")$x)
	if(!is.null(varying)){
		xName <- paste(xName, "X", deparse(attr(basis, "call")$varying), sep="")
		X <- cbind(varying, X * varying)
		Z <- Z * varying
	}	
	
	
	if(!is.null(by)){
		X.o <- X
		Z.o <- Z
		byName <- deparse(attr(basis, "call")$by)
		if(!allPen){
			basis$X <- basis$Z <- vector(mode="list", nlevels(by))
			for(i in 1:nlevels(by)){
				if(bySetToZero){
					keep <- 1*(by == levels(by)[i])
				} else keep <- rep(1, length(by))	
				#set X,Z partially to zero for each level of by-variable
				if(NCOL(X.o)){ 
					basis$X[[i]] <- X.o * keep
					#naming scheme: x.fx.bylevel1.fx1, x.fx.bylevel1.fx2, ...,x.fx.bylevel2.fx1, or xXvarying.fx.bylevel1.fx1
					
					colnames(basis$X[[i]]) <- paste(xName,".",byName,levels(by)[i],
							paste(".fx",1:NCOL(basis$X[[i]]),sep=""), sep="")
				}
				basis$Z[[i]] <- Z.o * keep
			}
			#naming scheme: bylevel1, bylevel2, ...
			names(basis$X) <- names(basis$Z) <- paste(byName, levels(by), sep="")
		} else {
			basis$X <- basis$Z <- vector(mode="list", 1)
			by <- C(by[, drop=TRUE], contr.treatment) #make sure treatment contrasts are used, unused levels dropped
			#basis$Z[[1]] <- Matrix(model.matrix(~ 0 + Z.o:by)) #TODO: can this be done without intermediate dense matrix?
			basis$Z[[1]] <- model.matrix(~ 0 + Z.o:by) #FIXME: ?constructing directly as sparse breaks transposing Z in subAZ?
			#cbind Z set partially to zero for each level of by-variable:
			## for(i in 1:nlevels(by[, drop=TRUE])) {
			##     if(bySetToZero){
			##         keep <- (by == levels(by[, drop=TRUE])[i])
			##     } else keep <- 1
			##     basis$Z[[1]] <- cBind(basis$Z[[1]], Z.o * keep)
			##     #basis$Z[[1]] <- Matrix(model.matrix(~ 0 + Z.o:by[, drop=TRUE]))
			## }	
			basis$X[[1]] <- X.o
			if(NCOL(X.o)) colnames(basis$X[[1]]) <- paste(xName,".",byName,paste(".fx",1:NCOL(X),sep=""), sep="")
			#naming scheme: u.x.by (also: name of duplicated by-variable in expandMf, subFcts)
			names(basis$X) <- paste("u", xName, byName, sep=".")
		}
	} else {
		if(NCOL(X)) colnames(X) <- paste(xName,paste(".fx",1:NCOL(X),sep=""), sep="")
		basis$X <- list(X)
		basis$Z <- list(Z)
	}  
	
	return(basis)
}

