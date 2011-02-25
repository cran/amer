safeDeparse <- function(expr){
	ret <- paste(deparse(expr), collapse="")
	#rm whitespace
	gsub("[[:space:]][[:space:]]+", " ", ret)
}

