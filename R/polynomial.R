constructConstPolyList=function(coefs)
{
	ans=do.call(polylist, as.list(coefs)) 
	class(ans)='mpolyList' ## class polylist
	ans
}
constructLinearPolyList=function(const.coefs=0, linear.coefs=rep(0, length(const.coefs)), linear.name='x')
{
	stopifnot(length(const.coefs) == length(linear.coefs))
	
	const = constructConstPolyList(const.coefs)
	idx = which(linear.coefs != 0)
	if(length(idx)==0L) return(const)
	
	oneX = polynomial(0:1)
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
	stopifnot(all(sapply(x, is.polynomial)))
	class(x)='mpolyList'
	x
}
as.mpolyList.polynomial = function(x, ...)
{
	structure(list(x), class='mpolyList')
}


'^.mpolyList' = function(e1, e2)
{
	stopifnot(is.numeric(e2) && length(e2)==1L)
	structure(lapply(e1, "^", e2), class="mpolyList")
}
	'*.mpolyList' = function(e1, e2) ## needed for polynom package
	{
		if(is.numeric(e2) || is.polynomial(e2)) {
			e2 = as.polynomial(e2)
			structure(lapply(e1, "*", e2), class="mpolyList")
		}else if(inherits(e2, 'mpolyList')){
			l1=length(e1); l2=length(e2)
			if(l1!=l2){
				L=max(c(l1,l2))
				e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
			}
			structure(mapply("*", e1, e2, SIMPLIFY=FALSE), class='mpolyList')
		}else stop('unknown e2 class')
	}
	'+.mpolyList' = function(e1, e2) ## needed for polynom package
	{
		if(is.numeric(e2) || is.polynomial(e2)) {
			e2 = as.polynomial(e2)
			structure(lapply(e1, "+", e2), class="mpolyList")
		}else if(inherits(e2, 'mpolyList')){
			l1=length(e1); l2=length(e2)
			if(l1!=l2){
				L=max(c(l1,l2))
				e1=rep(e1, length.out=L); e2=rep(e2, length.out=L)
			}
			structure(mapply("+", e1, e2, SIMPLIFY=FALSE), class='mpolyList')
		}else stop('unknown e2 class')
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
	mpl=c(...)
	if(length(mpl) == 1L) return(mpl[[1L]])
	Reduce("+", mpl[-1L], mpl[[1L]])
}



as.qmpolyList=function(x, ...) UseMethod("as.qmpolyList")
as.qmpolyList.mpolyList=function(x, denom.mpolyList=as.mpolyList(rep(1,length(x))), ...)
{
	structure(list(numerator=x, denominator=denom.mpolyList), class='qmpolyList')
}
commonDenom=function(numList, denomList, ...) ## need to optimize for speed
{  
	n=length(numList)
	stopifnot(n ==length(denomList))
	if(n==1) return(as.qmpolyList(numList, denomList))
	num.ans = polynomial(0)
	for(i in seq(n)){
		num.ans = num.ans + numList[[i]] * product(as.mpolyList(denomList[-i]))
	}
	denom.ans = product(denomList)
	as.qmpolyList(as.mpolyList(num.ans), as.mpolyList(denom.ans))
}

"+.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(0)))
{
	commonDenom(as.mpolyList(c(e1$numerator, e2$numerator)), 
				as.mpolyList(c(e1$denominator, e2$denominator)))
}
"-.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(0)))
{
	commonDenom(as.mpolyList(c(e1$numerator, e2$numerator * -1 )), 
				as.mpolyList(c(e1$denominator, e2$denominator)))
}
summation.qmpolyList=function(...)
{
	qmpl = list(...)
	ans0 = "+.qmpolyList"(qmpl[[1L]])
	if(length(qmpl)==1L) return(ans0)
	Reduce("+.qmpolyList", qmpl[-1L], ans0)
}
"*.qmpolyList"=function(e1, e2=as.qmpolyList(as.mpolyList(1)))
{
	e1=summation(e1); e2=summation(e2)
	as.qmpolyList(e1$numerator * e2$numerator, e1$denominator * e2$denominator)
}

degree=function(x,...)UseMethod('degree')
degree.polynomial=function(x,...)
{
	seq_along(x)-1L
}
#coef.polynomial is already defined in polynom
