trimZeros=function(vec)
{
  nzeros=sum(cumsum(abs(rev(vec)))==0)
  if(nzeros>0) vec[-seq(length(vec), length=nzeros, by=-1L)] else vec
}

eigenPolyRoot=function(co, only.real=FALSE)
{
  co=trimZeros(co)
  N=length(co)
  mat=Matrix::Matrix(0, N-1L, N-1L, sparse=TRUE)
  if(N==3){
	mat[2L,1L]=1
  }else if(N==2L){
	ans = -co[1L]/co[2L]
	return(if(isTRUE(only.real)) ans else as.complex(ans))
  }else if(N==1L){
	return(if(isTRUE(only.real)) numeric(0) else complex(0))
  }else diag(mat[-1L,])=1
  mat[,N-1L]=-co[-N]/co[N]
  rts=eigen(mat, symm=FALSE, only=TRUE)$value ## actually eigen converts mat to a dense matrix
  rts=rts[order(abs(Im(rts)), -(Re(rts)))]
  if(isTRUE(only.real)) Re(rts[Im(rts)==0]) else rts
}


numPolyVar.default=function(e1, ...)
{# numeric or bigq
	co=sign(e1[e1!=0])
	sum(co[-1]*co[-length(co)]<0)
}
numPolyVar.list=function(e1, at, ...)
{
	e1=do.call(c, lapply(e1, numPolyEval, at=at))
	NextMethod('numPolyVar')	
}
numPolyVar=function(e1, ...) UseMethod('numPolyVar')

numPolyRootUBound=function(e1, method=c('Cauchy', 'Lagrange', 'Kojima','Fujiwara','SumAdjRatio'))
{
	e1=trimZeros(e1)
	n=length(e1)
	method=match.arg(method, several.ok=TRUE)
	nmethod=length(method)
	bnd=vector('list', nmethod)
	names(bnd)=method
	if('Cauchy'%in% method)	bnd$Cauchy = 1 + max (abs(e1[-n]/e1[n]))
	if('Lagrange'%in% method) bnd$Lagrange = max(1, sum(abs(e1[-n]/e1[n])))
	if('Kojima'%in% method) bnd$Kojima = if(any(e1==0)) Inf else 2 * max(abs(e1[-n]/e1[-1L])*c(.5, rep(1,max(0,n-2))))
	if('Fujiwara'%in% method) bnd$Fujiwara = 2 * max((abs(e1[-n]/e1[n])*c(.5, rep(1,max(0,n-2))))^(1/safeseq(n-1,1,by=-1)))
	if('SumAdjRatio'%in% method) bnd$SumAdjRatio = if(any(e1[-1]==0)) Inf else sum(abs(e1[-n]/e1[-1]))
	ans=min(unlist(bnd), Inf, na.rm=TRUE)
	attr(ans, 'bounds')=bnd
	ans
}

numPolyRootLBound=function(e1, method=c('Rouche'))
{
	e1=trimZeros(e1)
	n=length(e1)
	method=match.arg(method, several.ok=TRUE)
	nmethod=length(method)
	bnd=vector('list', nmethod)
	names(bnd)=method
	if('Rouche'%in% method) bnd$Rouche = max(abs(e1[1L])/(abs(e1[1L]) + max(abs(e1[-1L]))), 
											 abs(e1[1L])/max(abs(e1[1L]),sum(abs(e1[-1L]))))
	ans=max(0, unlist(bnd), na.rm=TRUE)
	attr(ans, 'bounds')=bnd
	ans
}

numPolySturmSeq=function(e1)
{	if(is.numeric(e1[[1L]])) {
		neg1=-1 ; zero=0
	}else {
		neg1=as.bigq(-1); zero=as.bigq(0)
	}
	e1=trimZeros(e1)
	dif=numPolyDeriv(e1)
	ans=list(e1, dif)
	last=dif; last2=e1
	repeat{
		tmp = numPolyMult(neg1, numPolyDiv(last2, last)$remainder)
		if(length(tmp)==0) tmp=zero
		ans=c(ans, list(tmp))
		if(length(tmp)==1L) break
		#if(length(tmp)==0L) stop('not square free?')
		last2=last
		last=tmp
	}
	ans
}

numPolyNRealRoot=function(e1, lower=-upper, upper=numPolyRootUBound(e1)*1.001, method='Sturm')
{	stopifnot(upper>lower)
	sqFree=numPolySquareFree(e1)
	sqFree=sqFree[sapply(sqFree,length)>1L]
	
	numPolyNRealRoot.Recall=eval(match.call(), parent.frame())
	if(length(sqFree)>1L) return(sum(sapply(sqFree, numPolyNRealRoot.Recall, lower,upper, method)))
	
	e1=sqFree[[1L]]
	method=match.arg(method)
	if(method!='Sturm') .NotYetImplemented()
	### sturm's theorem
		sturm.seq=numPolySturmSeq(e1)
		return(numPolyVar(sturm.seq, lower) - numPolyVar(sturm.seq, upper))
}

numPolyRealRootIso=function(e1, lower=-upper, upper=numPolyRootUBound(e1)*1.001, eps=1e-3)
{## bisection based on Sturm's theorem
	if(is.infinite(lower)) lower = -numPolyRootUBound(e1)*1.001
	if(is.infinite(upper)) upper =  numPolyRootUBound(e1)*1.001

	nroots=numPolyNRealRoot(e1, lower, upper)
	if(nroots==0) return(list())

	if(upper - lower < eps && nroots==1) return(list(c(lower, upper)))
	
	mid=.5*(lower+upper)
#	ans=list()
#	if(numPolyNRealRoot(e1, lower, mid) >=1) ans=c(ans, Recall(e1, lower, mid, eps))
#	if(numPolyNRealRoot(e1, mid, upper) >=1) ans=c(ans, Recall(e1, mid, upper, eps))
#	return(ans)
	numPolyRealRootIso.Recall=eval(match.call()[[1]], parent.frame()) ## for profiling purposes
	c(numPolyRealRootIso.Recall(e1, lower, mid, eps),numPolyRealRootIso.Recall(e1, mid, upper, eps))
}
