addSymbol=as.symbol('+')
multSymbol=as.symbol('*')
subtSymbol=as.symbol('-')

symPolyMult=function(e1, e2)
{
	d1=length(e1); # l1=sapply(e1, length)
	d2=length(e2); # l2=sapply(e2, length)
	ans=rep(expression(0), d1+d2-1L)

	for(i1 in seq(d1)-1L){
		j1=i1+1L
		for(i2 in seq(d2)-1L){
			j2=i2+1L
			j=i1+i2+1L
			toAdd=as.call(list(multSymbol, e1[[j1]], e2[[j2]]))
			ans[[j]]=as.call(list(addSymbol, ans[[j]], toAdd))
		}
	}
	ans
}
#symPolyMult(expression(0,1), expression(2,3))
#symPolyMult(expression(0,1), expression(1+1,x^2))
#symPolyMult(expression(0,1, sqrt(y)), expression(1+1,x^2))


symDot2List=function(...)
{
  dot.list=list(...)
  lists=list();
  for(i in seq_along(dot.list)){
    if(is.list(dot.list[[i]]) && is.expression(dot.list[[i]][[1L]])) {
      lists=c(lists, dot.list[[i]])
    }else if(is.expression(dot.list[[i]])){
      lists=c(lists, list(dot.list[[i]]))
    }else stop('... needs to be expressions or lists of expresions')
  }
  lists
}

symPolyProd=function(...)
{ ## This function computes:  Reduce(symPolyMult, symDot2List(...))
	lst=symDot2List(...)
	m=length(lst)
	if(m==1L) return(lst[[1L]])
	if(m==2L) return(symPolyMult(lst[[1L]], lst[[2L]]))
	
	ans=rep(list(expression()), ceiling(m/2))
	for(i in safeseq(1, m, by=2)){
		ans[[(i+1)/2]]=if(i==m) lst[[m]] else symPolyMult(lst[[i]], lst[[i+1L]])
	}
	Recall(ans)
}

symPolyPow=function(e1, pow)
{
  stopifnot(pow==round(pow))
  symPolyProd(rep(list(e1),pow))
}
#sapply(symPolyPow(expression(0,1,2), 2),eval)

symPolyAdd=function(e1, e2)
{
	d1=length(e1); # l1=sapply(e1, length)
	d2=length(e2); # l2=sapply(e2, length)
	if(d1<d2){tmp=e1; e1=e2; e2=tmp; d2=d1}
    for(i in seq(d2)){
        e1[[i]]=as.call(list(addSymbol, e1[[i]], e2[[i]]))
    }
    e1
}
#symPolyAdd(expression(0,1), expression(2,3))
#symPolyAdd(expression(0,1), expression(1+1,x^2))
#symPolyAdd(expression(0,1, sqrt(y)), expression(1+1,x^2))


symPolySum=function(...)
{#  Reduce(symPolyAdd, symDot2List(...))
	lst=symDot2List(...)
	m=length(lst)
	if(m==1L) return(lst[[1L]])
	if(m==2L) return(symPolyAdd(lst[[1L]], lst[[2L]]))
	
	ans=rep(list(expression()), ceiling(m/2))
	for(i in safeseq(1, m, by=2)){
		ans[[(i+1)/2]]=if(i==m) lst[[m]] else symPolyAdd(lst[[i]], lst[[i+1L]])
	}
	Recall(ans)
}

symPolySubt=function(e1, e2)
{
    if(missing(e2)){
        stopifnot(length(e1)>0L)
        e2=e1
        e1=expression(0)
    }
	d1=length(e1); 
	d2=length(e2); 
    if(d1<d2){
        e1=c(e1, rep(expression(0), d2-d1))
    }else if(d1>d2)         e2=c(e2, rep(expression(0), d1-d2))
    for(i in seq_along(e1)){
        e1[[i]]=as.call(list(subtSymbol, e1[[i]], e2[[i]]))
    }
    e1
}
#symPolySubt(expression(0,1), expression(2,3))
#symPolySubt(expression(0,1), expression(1+1,x^2))
#symPolySubt(expression(0,1, sqrt(y)), expression(1+1,x^2))

symPolyRatSimSum=function(numer, denom)
{
  if(!is.list(numer)) numer=list(numer)
  if(!is.list(denom)) denom=list(denom)
  n=length(numer); stopifnot(length(denom)==n)
  if(n==1) return(list(numerator=numer[[1L]], denominator=denom[[1L]]))
  
  ans.numer=numer
  for(i in seq(n)){
    ans.numer[[i]]=symPolyMult(numer[[i]], symPolyProd(denom[-i]))
  }
  list(numerator=symPolySum(ans.numer), 
       denominator=symPolyProd(denom))
}
#symPolyRatSimSum(list(expression(0,1),expression(2)),rep(list(expression(3,2,1)),2))

