subAZ <-
function(m, fct)
# replace design matrices for fake factors with designs for penalized spline basis
{
	for(i in 1:length(fct)){
		if(length(fct[[i]]$Z) == 1){
			ind <- which(names(m$FL$fl)[attr(m$FL$fl, "assign")]==names(fct[i])) 
			Zt <- as(t(fct[[i]]$Z[[1]]), "sparseMatrix")
			m$FL$trms[[ind]]$A <- m$FL$trms[[ind]]$Zt <- Zt
			dimnames(m$FL$trms[[ind]]$ST) <- list(deparse(attr(fct[[i]], "call")[[1]]), deparse(attr(fct[[i]], "call")[[1]])) 
		} else {
			for(j in 1:length(fct[[i]]$Z)){
				ind <- grep(paste(names(fct[i]),names(fct[[i]]$Z)[j],sep='.'), names(m$FL$fl)[attr(m$FL$fl, "assign")])
				Zt <- as(t(fct[[i]]$Z[[j]]), "sparseMatrix")
				m$FL$trms[[ind]]$A <- m$FL$trms[[ind]]$Zt <- Zt
				dimnames(m$FL$trms[[ind]]$ST) <- list(deparse(attr(fct[[i]], "call")[[1]]), deparse(attr(fct[[i]], "call")[[1]]))
			}
		}
	}
	return(m)
}

