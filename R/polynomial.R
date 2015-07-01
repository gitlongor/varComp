constructConstPolyList=function(coefs)
{
	tmp=lapply(coefs, list)
	ans = lapply(lapply(tmp, function(zzz){names(zzz[[1]])='coef';zzz}),mpoly)
	class(ans)='mpolyList'
	ans
}
constructLinearPolyList=function(const.coefs=0, linear.coefs=rep(0, length(const.coefs)), linear.name='x')
{
	stopifnot(length(const.coefs) == length(linear.coefs))
	
	const = constructConstPolyList(const.coefs)
	idx = which(linear.coefs != 0)
	if(length(idx)==0L) return(const)
	
	oneX = mpoly(list(structure(rep(1,2L), names=c(linear.name, 'coef'))))
	lin = lapply( constructConstPolyList(linear.coefs), "*", oneX)
	ans = mapply("+", const, lin, SIMPLIFY=FALSE)
	class(ans)='mpolyList'
	ans
}

as.mpolyList = function(x, ...) UseMethod("as.mpolyList")
as.mpolyList.numeric=function(x, ...)
{
	do.call('constructLinearPolyList', list(const.coefs=x, ...))
}
as.mpolyList.list=function(x, ...)
{
	stopifnot(all(sapply(x, is.mpoly)))
	class(x)='mpolyList'
	x
}
as.mpolyList.mpoly = function(x, ...)
{
	structure(list(x), class='mpolyList')
}


'^.mpolyList' = function(e1, e2)
{
	structure(lapply(e1, "^", e2), class="mpolyList")
}

product=function(...)UseMethod("product")
product.default=base::prod
product.mpolyList = function(...)
{
	mpl=c(...)
	if(length(mpl) == 1L) return(mpl[[1L]])
	Reduce("*", mpl[-1L], mpl[[1L]])
}

summation=function(...)UseMethod("summation")
summation.default=base::sum
summation.mpolyList = function(...)
{
	mpl=list(...)
	if(length(mpl) == 1L) return(mpl[[1L]])
	Reduce("+", mpl[-1L], mpl[[1L]])
}


as.qmpolyList=function(x, ...) UseMethod("as.qmpolyList")
as.qmpolyList.mpolyList=function(x, denom.mpolyList=as.mpolyList(rep(1,length(x))), ...)
{
	structure(list(numerator=x, denominator=denom.mpolyList), class='qmpolyList')
}

"+.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(0)))
{
	commonDenom(as.mpolyList(c(e1$numerator, e2$numerator)), 
				as.mpolyList(c(e1$denominator, e2$denominator)))
}
"-.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(0)))
{
	commonDenom(as.mpolyList(c(e1$numerator, e2$numerator * as.mpolyList.mpoly(mp('-1')) )), 
				as.mpolyList(c(e1$denominator, e2$denominator)))
}
"*.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(1)))
{
	e1=summation(e1); e2=summation(e2)
	as.qmpolyList(e1$numerator * e2$numerator, e1$denominator * e2$denominator)
}
summation.qmpolyList=function(...)
{
	qmpl = list(...)
	ans0 = "+.qmpolyList"(qmpl[[1L]])
	if(length(qmpl)==1L) return(ans0)
	Reduce("+.qmpolyList", qmpl[-1L], ans0)
}

commonDenom=function(numList, denomList, ...) ## need to optimize for speed
{  
	n=length(numList)
	stopifnot(n ==length(denomList))
	if(n==1) return(as.qmpolyList(numList, denomList))
	num.ans = mpoly(list(c(coef=0)))
	for(i in seq(n)){
		num.ans = num.ans + numList[[i]] * product(as.mpolyList(denomList[-i]))
	}
	denom.ans = product(denomList)
	as.qmpolyList(as.mpolyList(num.ans), as.mpolyList(denom.ans))
}


degree=function(x,...)UseMethod('degree')
degree.mpoly=function(x,...)
{
	nms=setdiff(unlist(lapply(x,names)), 'coef')
	if(length(nms)==0L) {
		rep(0, length(x))
	}else if(length(nms)==1L){
		ans = sapply(x, '[',nms)
		ans[is.na(ans)]=0
		ans
	}else .NotYetImplemented()
}

coef.mpoly = function(x, ...)
{
	sapply(x, '[','coef')
}	
