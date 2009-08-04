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


