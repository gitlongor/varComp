minque <-
function(y, varcov, start=rep(0, length(varcov)), lower.bound=-Inf, restricted=TRUE)
{ 
	k=varcov; varRatio=start
  if(!isTRUE(restricted)) .NotYetImplemented()
  stopifnot(is.numeric(y) && is.list(k) && is.numeric(varRatio))
  nK=length(k)
  varRatio=c(rep(varRatio, length=nK),1)
  k[[nK+1L]]=diag(1,nrow(k[[1L]]))  
  
  LI=t(backsolve(chol(Reduce(`+`,mapply(`*`,k,varRatio,SIMPLIFY=FALSE))), k[[nK+1]]))
  LIk=lapply(lapply(k, `%*%`, LI), crossprod, LI)
  LIy=LI%*%y

  S=outer(LIk,LIk,function(a,b)sapply(mapply(`*`,a,b,SIMPLIFY=FALSE),sum))
  u=sapply(lapply(LIk,'%*%',LIy),crossprod,LIy)

  if(lower.bound<0 && is.infinite(lower.bound)){
    ans0=ginv(S)%*%u  #solve(S, u) ## FIXME: add error handler
  }else{
    # require(quadprog)
	Dmat=crossprod(S); dvec=crossprod(S,u); Amat=rbind(diag(1,nK),0);bvec=rep(lower.bound,nK)
	repeat{
		qp.rslt=try(solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=bvec, meq=0L, factorized=FALSE), silent=TRUE)
		if(!inherits(qp.rslt, 'try-error')) break
		if(attr(qp.rslt, 'condition')$message==
			'matrix D in quadratic function is not positive definite!'){
				Dmat=as.matrix(nearPD(crossprod(S))$mat)
		}else if(attr(qp.rslt, 'condition')$message==
			'constraints are inconsistent, no solution!'){ 
				#truncate unconstrained solution and rescale
				unconstr=ginv(S)%*%u 
				s.unconstr=sum(unconstr)
				unconstr[unconstr<lower.bound]=lower.bound
				qp.rslt=list(solution=unconstr / sum(unconstr) * s.unconstr)
				break
		}else break
	}
    ans0=qp.rslt$solution
  }
  ans=ans0[-nK-1L]/ans0[nK+1L]
  ans
}
