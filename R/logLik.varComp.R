logLik.varComp=function(object, FUNCTION=FALSE, ...)
{
	if(isTRUE(FUNCTION)){
		return(object$PREML.fun)
	}
	REML=object$control$REML
	if(!isTRUE(REML)) .NotYetImplemented()
	val=object$PREML
	vc=coef(object, "varComp")
	p=sum(vc>object$control$boundary.eps)
	
    N = attr(val, "nall") <- nobs(object)
    attr(val, "nobs") <- N - REML * length(coef(object, 'fixed'))
    attr(val, "df") <- p #+ length(coef(object[["modelStruct"]])) + 1
    class(val) <- "logLik"	
	val
}

