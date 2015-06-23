trimZeros=function(vec, eps=.Machine$double.eps^.5)
{
  #nzeros=sum(abs(cumsum(abs(rev(vec))))<=eps)
  nzeros=sum(cumsum(abs(rev.default(vec))>eps)==0)
  if(nzeros>0) vec[-seq(length(vec), length=nzeros, by=-1L)] else vec
}

numPolyMult=function(e1, e2)
{
	d1=length(e1); # l1=sapply(e1, length)
	d2=length(e2); # l2=sapply(e2, length)
	if(d1==0) return(e1)
	if(d2==0) return(e2)
	ans=rep(if(is.numeric(e1[[1L]])) 0 else as.bigq(0), d1+d2-1L)

	for(i1 in (1:d1)-1L){
		j1=i1+1L
		i2 = (1:d2)-1L

		j2 = i2+1L
		j =i1+i2+1L
		ans[j]=ans[j]+e1[j1]*e2[j2]
		if(FALSE){
			for(i2 in seq(d2)-1L){
				j2=i2+1L
				j=i1+i2+1L
				ans[[j]]=ans[[j]]+e1[[j1]]*e2[[j2]] 
			}
		}
	}
	ans
}
#numPolyMult(c(0,1), c(2,3))
#numPolyMult(c(0,1), c(1+1,3^2))
#numPolyMult(c(0,1, sqrt(5)), c(1+1,5^2))


numDot2List=function(...)
{
  dot.list=list(...)
  lists=list();
  for(i in seq_along(dot.list)){
    if(is.list(dot.list[[i]]) && 
        (is.numeric(dot.list[[i]][[1L]]) || class(dot.list[[i]][[1L]])=='bigq')
        ) {
      lists=c(lists, dot.list[[i]])
    }else if(is.numeric(dot.list[[i]])){
      lists=c(lists, list(dot.list[[i]]))
    }else stop('... needs to be numeric, bigq or lists theirof')
  }
  lists
}


numPolyProd=function(...)
{ ## This function computes:  Reduce(numPolyMult, numDot2List(...))
	lst=numDot2List(...)
	m=length(lst)
	if(m==1L) return(lst[[1L]])
	if(m==2L) return(numPolyMult(lst[[1L]], lst[[2L]]))
	
	ans=rep(list(c()), ceiling(m/2))
	for(i in safeseq(1, m, by=2)){
		ans[[(i+1)/2]]=if(i==m) lst[[m]] else numPolyMult(lst[[i]], lst[[i+1L]])
	}
	numPolyProd.Recall=eval(match.call()[[1]], parent.frame())
	numPolyProd.Recall(ans)
}

numPolyPow=function(e1, pow)
{
  stopifnot(pow==round(pow))
  numPolyProd(rep(list(e1),pow))
}
#sapply(numPolyPow(expression(0,1,2), 2),eval)

numPolyAdd=function(e1, e2)
{
	d1=length(e1); # l1=sapply(e1, length)
	d2=length(e2); # l2=sapply(e2, length)
	if(d1==0)return(e2)
	if(d2==0)return(e1)
	if(d1<d2){tmp=e1; e1=e2; e2=tmp; d2=d1}
	idx=safeseq(1L,d2,by=1L)
	e1[idx]=e1[idx]+e2
	if(FALSE){
		for(i in seq(d2)){
			e1[[i]]=e1[[i]] + e2[[i]]
		}
	}
    e1
}
#numPolyAdd(expression(0,1), expression(2,3))
#numPolyAdd(expression(0,1), expression(1+1,x^2))
#numPolyAdd(expression(0,1, sqrt(y)), expression(1+1,x^2))


numPolySum=function(...)
{#  Reduce(numPolyAdd, numDot2List(...))
	lst=numDot2List(...)
	m=length(lst)
	if(m==1L) return(lst[[1L]])
	if(m==2L) return(numPolyAdd(lst[[1L]], lst[[2L]]))
	
	ans=rep(list(expression()), ceiling(m/2))
	for(i in safeseq(1, m, by=2)){
		ans[[(i+1)/2]]=if(i==m) lst[[m]] else numPolyAdd(lst[[i]], lst[[i+1L]])
	}
	numPolySum.Recall=eval(match.call()[[1]], parent.frame())
	numPolySum.Recall(ans)
}

numPolySubt=function(e1, e2)
{	zero=if(is.numeric(e1[[1L]])) 0 else as.bigq(0)
    if(missing(e2)){
        stopifnot(length(e1)>0L)
        e2=e1
        e1=zero
    }
	d1=length(e1); 
	d2=length(e2); 
    if(d1<d2){
        e1=c(e1, rep(zero, d2-d1))
    }else if(d1>d2)         e2=c(e2, rep(zero, d1-d2))

	if(FALSE){
		for(i in seq_along(e1)){
			e1[[i]]= e1[[i]] - e2[[i]]
		}
		e1
	}
	e1-e2
}
#numPolySubt(expression(0,1), expression(2,3))
#numPolySubt(expression(0,1), expression(1+1,x^2))
#numPolySubt(expression(0,1, sqrt(y)), expression(1+1,x^2))

numPolyDiv=function(e1, e2)
{	zero=if(is.numeric(e1[[1L]])) 0 else as.bigq(0)
	e1=trimZeros(e1)
	e2=trimZeros(e2)
	
	q=zero; r=e1
	d=length(e2)-1; C=tail(e2,1L)
	deg.r=length(r)-1
	while (deg.r >= d){
		s = c(rep(zero, deg.r-d), tail(r,1L)/C)
		q = trimZeros( numPolyAdd(q, s) )
		r = trimZeros( numPolySubt(r, numPolyMult(s, e2)) )
		deg.r = length(r) -1
	}
	list(quotient = q, remainder = r)
}

numPolyGCD=function(e1, e2, trace=FALSE)
{
	if(isTRUE(trace)) print(e2)
	if(length(trimZeros(e2))==0L) return(e1)
	numPolyGCD.Recall=eval(match.call()[[1]], parent.frame())
	ans=(numPolyGCD.Recall(e2, numPolyDiv(e1, e2)$remainder, trace))
	ans/ans[length(ans)]  ## monic
}

numPolySquareFree=function(e1)
{
	e1=trimZeros(e1)
	#if(length(e1)==1L) return(e1/e1)
	#e1=e1/e1[length(e1)]
	
	a0=numPolyGCD(e1, numPolyDeriv(e1))
	b1=numPolyDiv(e1, a0)$quotient
	c1=numPolyDiv(numPolyDeriv(e1), a0)$quotient
	d1=numPolySubt(c1, numPolyDeriv(b1))
	i=1
	ans=list()
	repeat{
		ai=numPolyGCD(b1, d1)
		b2=numPolyDiv(b1, ai)$quotient
		c2=numPolyDiv(d1, ai)$quotient
		ans=c(ans, list(ai))
		i=i+1
		d1=numPolySubt(c2, numPolyDeriv(b2))
		b1=trimZeros(b2); c1=c2
		if(length(b1)==1L && b1==1) break
	}
	ans
}
