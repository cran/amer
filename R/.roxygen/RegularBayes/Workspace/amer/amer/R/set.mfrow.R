set.mfrow <-
function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE)
# taken from coda 
{
	mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
				if (Nchains == 2) {
					switch(min(Nparms, 5), c(1, 2), c(2, 2), c(3, 2), 
							c(4, 2), c(3, 2))
				}
				else if (Nchains == 3) {
					switch(min(Nparms, 5), c(2, 2), c(2, 3), c(3, 3), 
							c(2, 3), c(3, 3))
				}
				else if (Nchains == 4) {
					if (Nparms == 1) 
						c(2, 2)
					else c(4, 2)
				}
				else if (any(Nchains == c(5, 6, 10, 11, 12))) 
					c(3, 2)
				else if (any(Nchains == c(7, 8, 9)) || Nchains >= 13) 
					c(3, 3)
			}
			else {
				if (nplots == 1) {
					mfrow <- switch(min(Nparms, 13), c(1, 1), c(1, 2), 
							c(2, 2), c(2, 2), c(3, 2), c(3, 2), c(3, 3), 
							c(3, 3), c(3, 3), c(3, 2), c(3, 2), c(3, 2), 
							c(3, 3))
				}
				else {
					mfrow <- switch(min(Nparms, 13), c(1, 2), c(2, 2), 
							c(3, 2), c(4, 2), c(3, 2), c(3, 2), c(4, 2), 
							c(4, 2), c(4, 2), c(3, 2), c(3, 2), c(3, 2), 
							c(4, 2))
				}
			}
	return(mfrow)
}

