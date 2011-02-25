expandMf <-
function(fr, fct)
# cbind model frame with design matrices for unpenalized&penalized parts of the smooth fcts.  
{
	for(i in 1:length(fct)){
		#matrix with all unpenalized terms for fct
		newX <- do.call(cBind, fct[[i]]$X)
		
		#factor variables with no. of levels = no. of penalized basis fcts  
		#newFact <-   replicate(length(fct[[i]]$Z), rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]])))
		newFact <- data.frame(factor(rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]]))))
		if(length(fct[[i]]$Z) > 1){
			for(j in 2:length(fct[[i]]$Z)){
				newFact <- cbind(newFact, factor(rep(1:ncol(fct[[i]]$Z[[1]]), length=nrow(fct[[i]]$Z[[1]]))))
			}
		}
		
		colnames(newFact) <- if(length(fct[[i]]$Z) == 1){ 
					names(fct)[i]
				} else {
					paste(names(fct)[i],".",names(fct[[i]]$Z),sep="")
				}
		
		if(eval(attr(fct[[i]],"call")$allPen)){
			# duplicate grouping factor for allPen-function groups so that assignment (which entries in ranef belong to which penalization 
			# group) can be reconstructed from the fitted model object m if there is another random effect associated with the by-variable.
			# will need this for predict etc.. since attr(m@flist,"assign") only works the other way around....
			newFact <- cBind(newFact, eval(attr(fct[[i]],"call")$by, fr))
			colnames(newFact)[ncol(newFact)] <- names(fct[[i]]$X)
		} 
		
		fr <- cBind(cBind(fr, newX),newFact)
	}
	return(fr)
}

