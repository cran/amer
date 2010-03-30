subFcts <-
function(rhs, fctterm, fct, fr)
# replace formula parts for smooth functions with  xi + (xi^2+ )... + xi^dimUnpen + (1|fcti) or
# by*(xi + xi^2+ ... + xi^dimUnpen) + (1|fcti.1) + ... + (1|fcti.N) for by-variable with N levels	
{
	for(i in 1:length(fct)){
		by <- eval(attr(fct[[i]],"call")$by, fr)
		allPen <- eval(attr(fct[[i]],"call")$allPen)
		diag <- eval(attr(fct[[i]],"call")$diag)

		replacement <-  
				if(is.null(by)){
					# 1 + x.fx1 + x.fx2+ ... + (1|f.x)
					paste(ifelse(ncol(fct[[i]]$X[[1]])!=0,
								paste(as.vector(sapply(fct[[i]]$X,colnames)),collapse=" + "),
								"1"),
							" + (1|",names(fct)[i],")",sep="")		
				} else {
					if(allPen){ 
						if(!diag){
							# add correlated random effects for normally unpenalized part of basis grouped according to by and fake random intercept
							# (1 + x.fx1 + x.fx2+ ...|u.x.by)  + (1|f.x.by)
							paste(
							 paste(paste("(1",
										paste(as.vector(sapply(fct[[i]]$X,colnames)), collapse="+"),
									sep="+"),
							  "|", 
							  names(fct[[i]]$X),")",
							  sep=""),
							  	paste("(1|",names(fct)[i],")",sep="", collapse=" + "),
								sep =" + ")
						} else {
							# add independent random effects for normally unpenalized part of basis grouped according to by and fake random intercept
							# (1|u.x.by) + x.fx1|u.x.by) + x.fx2|u.x.by) + ...  + (1|f.x.by)
							paste(
									paste(c("(1", paste("(0+", as.vector(sapply(fct[[i]]$X,colnames)),sep="")),"|", names(fct[[i]]$X),")",sep="",collapse=" + "),
									paste("(1|",names(fct)[i],")",sep="", collapse=" + "),
									sep =" + ")
						}
								
					} else { 
						# add fixed effect for unpenalized part of basis + fake random intercept for each by-level
						# by + x.fx1.BYlevel1 + x.fx2.BYlevel1 +...+ (1|f.x.BYlevel1) + ... + x.fx1.BYlevelD + x.fx2.BYlevelD +... + (1|f.x.BYlevelD)
						paste(#deparse(attr(fct[[i]],"call")$by), 
								ifelse(ncol(fct[[i]]$X[[1]])!=0,
										paste(as.vector(sapply(fct[[i]]$X,colnames)),collapse=" + "),
										"1"),
								paste("(1|",names(fct)[i],".",names(fct[[i]]$Z),")",sep="", collapse=" + "),
								sep =" + ")
					}
				}
		rhs <- sub(deparse(fctterm[[i]]), replacement, rhs, fixed=T)
	}
	return(rhs)
}

