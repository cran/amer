indsF <-
function(m, fct, fctterm){
# add assign-like info to fctterm:
# which penalization/ranef groups and coefficients (fixed/random) belong to which function 
# also include info on global intercept and by-level intercepts	
	ranefinds <- lme4:::reinds(m@Gp)
	
	indIntercept <- ifelse("(Intercept)" %in% names(fixef(m)), 1, 0)
	
	for(i in 1:length(fctterm)){
		if(length(fct[[i]]$Z) == 1){
			
			attr(fctterm[[i]], "indGrp") <- match(names(fct)[i], colnames(m@flist)) 
			if(eval(attr(fct[[i]], "call")$allPen)) {
				#add pen. group(s) with grouping factor u.x.by
				indUGrp <- match(sub("f.", "u.", names(fct)[i]), colnames(m@flist))
				attr(fctterm[[i]], "indGrp") <-  c(attr(fctterm[[i]], "indGrp"), which(attr(m@flist, "assign")==indUGrp))
			}
			attr(fctterm[[i]], "indPen") <- unlist(ranefinds[attr(fctterm[[i]], "indGrp")])
			
			if(!(eval(attr(fct[[i]], "call")$allPen)||ncol(fct[[i]]$X[[1]])==0)){
				attr(fctterm[[i]], "indUnpen") <-  sapply(paste("^",colnames(fct[[i]]$X[[1]]),"$",sep=""),
						grep, x=names(m@fixef))
				names(attr(fctterm[[i]], "indUnpen")) <- colnames(fct[[i]]$X[[1]])
			} else attr(fctterm[[i]], "indUnpen") <- 0

			attr(fctterm[[i]], "indConst") <- indIntercept
					
			attr(fctterm[[i]], "indGrp") <- list(attr(fctterm[[i]], "indGrp"))
			attr(fctterm[[i]], "indPen") <- list(attr(fctterm[[i]], "indPen"))
			attr(fctterm[[i]], "indUnpen") <- list(attr(fctterm[[i]], "indUnpen"))
			attr(fctterm[[i]], "indConst") <- list(attr(fctterm[[i]], "indConst"))
		} else {
			by <- eval(attr(fct[[i]],"call")$by, m@frame)
			attr(fctterm[[i]], "indGrp") <- vector(mode="list", length=nlevels(by))
			attr(fctterm[[i]], "indPen") <-	vector(mode="list", length=nlevels(by))
			attr(fctterm[[i]], "indUnpen") <- vector(mode="list", length=nlevels(by))
			attr(fctterm[[i]], "indConst") <- vector(mode="list", length=nlevels(by))
			for(j in 1:nlevels(by)){
				attr(fctterm[[i]], "indGrp")[[j]] <- grep(paste("^",paste(names(fct)[i],".",names(fct[[i]]$Z)[j],sep=""), "$", sep=""), colnames(m@flist))
				attr(fctterm[[i]], "indPen")[[j]] <- ranefinds[[attr(fctterm[[i]], "indGrp")[[j]]]]
				if( ncol(fct[[i]]$X[[j]]) == 0){
					attr(fctterm[[i]], "indUnpen")[[j]] <-	 0
				} else {
					attr(fctterm[[i]], "indUnpen")[[j]] <- sapply(
							paste("^",colnames(fct[[i]]$X[[j]]),"$",sep=""), grep, x=names(m@fixef))
					names(attr(fctterm[[i]], "indUnpen")[[j]]) <- colnames(fct[[i]]$X[[j]])
				}	
				#add by-level intercept:
				indBy <- grep(paste("^",safeDeparse(attr(fct[[i]],"call")$by), levels(by)[j],"$", sep=""), names(m@fixef))
				indBy <- indBy[!(indBy %in% attr(fctterm[[i]], "indUnpen")[[j]])]
				attr(fctterm[[i]], "indConst")[[j]] <- c(indIntercept, indBy)
			}
		}
	}
	return(fctterm)
}

