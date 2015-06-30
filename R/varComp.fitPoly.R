varComp.fitPoly = function(Y, X=matrix(0,length(Y),0L), K, control=varComp.control())
{
	{
	#                  a.	varComp: The main fitting function using REML criterion
	#                      i.	Y: response
	#                      ii.	X: fixed effect design matrix. Treated as intercept if missing. A zero matrix can be supplied for REML. 
	#                      iii.	Z: random effect design matrix, or a list of such matrices. Treated as identity if missing. 
	#                      iv.	K: A matrix or a list of matrices, each determining the correlations among observations. This has to be at least p.s.d. in this implementation. 
	#                      v.	information: A character specifying the information matrix to be used during optimization. For doing score tests, this has to be EI; for fitting models, others are OK and the default AOI is usually fast enough.
	#                      vi.	optMethod: A character specifying the optimization method. The default should work well in most cases. 'NRGD' is only a naive implementation and has not been tested extensively. 
	#                      vii.	boundary.eps: A small numeric, below which estimates of tau will be checked for hitting boundary (0). This is a simple attempt at differentiating small but non-zero variance components vs truly zero variance components. Another way of doing this is to divide the K matrices by some large number. 
	#                      viii.	conv: A small numeric, the convergence criterion, only used when optMethod='NRGD'.
	#                      ix.	nStepHalving: A positive integer, giving the max number of step halving during NRGD. 
	#                      x.	nIter: A positive integer, giving the max number of iterations allowed for NRGD.
	#                      xi.	starts: A vector of nonnegative doubles of the same length as the number of variance components. This is the starting value for tau, i.e., the ratio of each variance components to the error variance. 
	#                      xii.	plot.it: A logical scalar, indicating if PREML surface will be plotted (for one variance components only). Please send Long Qu an Email with your data set if you find multiple local maxima, for possible improvements of this function. 
	#                      xiii.	verbose: A logical scalar, indicating if extra information is printed. 
	#                             keepXYK: Logical, indication if original X, Y and K matrices are stored in the results. 
	}
  #if(!is.matrix(X)) X=as.matrix(X)  ## does it matter?
  #if(diff(range(X))>0) X=model.matrix(~X)  ##  dealing with constant X including zero X
  #Y=as.numeric(as.vector(Y)) ## plan to move to varComp formula interface for efficiency
	{
		 optMethod=control$optMethod
		 verbose= control$verbose 
		 starts= control$starts 
		 REML = control$REML 
			if(!isTRUE(REML)) stop("Currently only REML method is implemented")
		 information = control$information 
		 boundary.eps= control$boundary.eps 
		 nlminb.control = control$nlminb 
		 plot.it= control$plot.it 
		 keepXYK= control$keepXYK
		 
		 environment(V)=
		 environment(updateLI)=
		 environment(updateLIkLI)=
		 environment(updateLIy)=
		 environment(PREML)=
		 environment(preprocPREML)=
		 environment(obj)=
		 environment(obj2)=
		 environment(updateNumsPart)=
		 environment(updateNums)=
		 environment(updateDenom)=
		 environment(updateTr1)=
		 environment(updateTr2)=
		 environment(updateNums2)=
		 environment(score)=
		 environment(updateNegGrad)=
		 environment(gradFunc)=
		 environment(OI)=
		 environment(hess)=
		 environment(EI)=
		 environment(AOI)=
		 environment(AEI)=
		 environment(WAI)=
		 environment(updateNegHess)=sys.frame(sys.nframe())
	}

  
  qrx=qr(X)
  # * this block checks if an intercept is in the model or not; 
  # * don't needed this anymore
	 # Q1=qr.Q(qrx)
	 # if(qrx$rank>0 && max(abs(Q1%*%crossprod(Q1, rep(1,length(Y))) - rep(1,length(Y)))) > sqrt(.Machine$double.eps)){  ## no intercept
	   # X=cbind(`(Intercept)`=1, X)
	   # qrx=qr(X)
	 # }
  
  Q2=if(qrx$rank>0) qr.Q(qrx, complete=TRUE)[, -seq_len(qrx$rank),drop=FALSE] else diag(1, length(Y)) #  qr.Q(qrx, complete=TRUE) == IdMat
  y=crossprod(Q2, Y)
  n=length(y)
  
  if(missing(K) || length(K)==0) {  ## linear model fit
		null.sig2=drop(crossprod(y))  / n   
		null.preml=.5*(
		  -n*log(crossprod(y)) -n -n*log(2*pi)-n*log(n)
		)

	  ans=list(
	    ## varComp.fit specific block
		parms=numeric(0L),
		gradients=numeric(0L), 
		hessian=matrix(numeric(0L), 0L, 0L),
		sigma2=drop(null.sig2), 
		varComps=numeric(0),
		n.iter=0L, PREML=drop(null.preml),
		X.Q2=Q2, 
		residual.contrast=y, 
		working.cor=vector('list', 0L),
		
		## varComp common block
		# na.action=NULL,
		# offset = NULL,
		# contrasts=NULL,
		# xzlevels = NULL,
		# terms = NULL, 
		call=match.call(), 
		nobs = length(Y), 
		control=control, 
		random.labels = character(0L), 
		doFit= TRUE,
		
		# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
		X=if(isTRUE(keepXYK)) X else NULL,
		# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
		Y=if(isTRUE(keepXYK)) Y else NULL, 
		K=if(isTRUE(keepXYK)) vector('list', 0L) else NULL
	  )
	  class(ans)='varComp'	
	  return(ans)
  }
  
  nK=length(K)
  if(nK>1) plot.it=FALSE
  
  k=vector('list',nK)
  for(j in 1:nK) k[[j]]=crossprod(Q2, K[[j]]%*%Q2)
  
  LIkLI=k; LIy=numeric(n); numsPart=matrix(NA_real_, n, nK); 
  tr2=negHess=nums2=matrix(NA_real_, nK, nK)
  tr1=nums=negGrad=numeric(nK); denom=numeric(1L)
  LI=matrix(NA_real_, n, n)
  diag.1.n=diag(1,n)
  if(nK==1) {
	eigK=eigen(k[[1]],TRUE)
	eigK$tvector=t(eigK$vector)
  }
  
  stdZ = drop(crossprod(eigK$vector , y))
  sum2.denom=vector('list',n)
  for(i in seq(n)){
	termList=list(c(coef=1),c(tau=1, coef=eigK$value[i]))
	sum2.denom[[i]]=mpoly(termList)
  }
  class(sum2.denom)='mpolyList'
  sum1.denom = lapply(sum2.denom,"^",2)
  sum3.denom = sum2.denom
  
  constructConstPolyList=function(coefs)
  {
	tmp=lapply(coefs, list)
    ans = lapply(lapply(tmp, function(zzz){names(zzz[[1]])='coef';zzz}),mpoly)
	class(ans)='mpolyList'
	ans
  }
  sum1.num = constructConstPolyList(eigK$value * stdZ^2 * n)
  sum2.num = constructConstPolyList(eigK$value)
  sum3.num = constructConstPolyList( stdZ^2 )
  
  prod.mpolyList = function(mpl)
  {
	if(length(mpl) == 1L) return(mpl)
	Reduce("*", mpl[-1L], mpl[[1L]])
  }
  
  commonDenom=function(numList, denomList)
  {
	n=length(numList)
	stopifnot(n ==length(denomList))
	num.ans = mpoly(list(c(coef=0)))
	for(i in seq(n)){
		num.ans = num.ans + numList[[i]] * prod.mpolyList(denomList[-i])
	}
	denom.ans = prod.mpolyList(denomList)
	list(numerator = num.ans, denominator=denom.ans)
  }
  sum2 = commonDenom(sum2.num, sum2.denom)
  sum3 = commonDenom(sum3.num, sum3.denom)
  sum1 = commonDenom(sum1.num, sum1.denom)
  sum23 = list(numerator = sum2$numerator * sum3$numerator, denominator = sum2$denominator * sum3$denominator)
  ans.num = sum1$numerator * sum23$denominator - sum23$numerator * sum1$denominator
  degree = sapply(ans.num, '[','tau'); degree[is.na(degree)]=0
  polyCoefs = sapply(ans.num, '[','coef')[order(degree)]

  roots = polyroot(polyCoefs)
  candidates = Re(roots)[Re(roots)>=0 & abs(Im(roots)) < .Machine$double.eps^.5]
  n.nr=0L
  
  infoFunc=get(information)
  preprocScore=expression({
	updateLIkLI()
	updateNumsPart();   
	updateNums(); updateDenom(); updateTr1()
  })
  preprocInfo=switch(information,
	OI=expression({ updateNums2(); updateTr2() }),
	EI=expression({ updateTr2();  }),
	AOI=expression({updateNums2();    }),
	AEI=expression({updateNums2();   }),
	WAI=expression({updateNums2();   })
  )
  
if(FALSE){
  last.tau=tau=if(is.null(starts)) {
	minque(y, k, rep(0,nK), lower.bound=.Machine$double.eps^.5, restricted=TRUE)
  }else rep(starts, length=nK)
  ltau=log(max(.Machine$double.eps, tau))
  
  if(optMethod=='nlminb'){
	objNeg=function(tau)-obj(tau)
	gradNeg=function(tau)-gradFunc(tau)
	hessNeg=function(tau)-hess(tau)
	nlminb.fit=nlminb(tau, objNeg, gradNeg, hessNeg, lower=rep(0,nK), control=nlminb.control)
	tau=nlminb.fit$par
	n.nr=nlminb.fit$iterations
  }else if(optMethod=='optim'){
	optim.fit=optim(tau, obj, gradFunc, method='L-BFGS-B', lower=rep(0, nK), control=list(fnscale=-1))
	tau=optim.fit$par
	n.nr=optim.fit$counts
  }else if(optMethod=='NRGD'){
	nStepHalving = 20L
	preprocPREML(tau)
	last.func = PREML()
	n.nr=0L
	repeat{
	  if(verbose) cat('PREML=',last.func, 'tau=',tau,'\n')
	  eval(preprocScore)
	  eval(preprocInfo)
	  updateNegGrad()
	  updateNegHess()
	  
	  negGrad=tau*negGrad
	  negHess=negHess*outer(tau,tau)
	  diag(negHess)=diag(negHess)+negGrad
	  
	  adj.ltau=solve(negHess, negGrad)
	  
	  doNewton=TRUE
	  repeat{
		step.size=1
		n.step=0L
		ltau.new=ltau
		repeat{
		  ltau.new=ltau - step.size * adj.ltau
		  preprocPREML(exp(ltau.new))
		  this.func=PREML()
		  if( ( (this.func>last.func) & doNewton) | ((!doNewton) & (this.func>=last.func)) ) {
			doNewton=TRUE
			break
		  }
		  step.size=step.size/2
		  n.step=n.step+1L
		  if((n.step > nStepHalving) & doNewton ) {  ## try gradient descent
			#warning('Step halving failed')
			adj.ltau=negGrad
			doNewton=!doNewton
			break
		  }
		}
		if(isTRUE(doNewton)) break
	  }
	  
	  last.tau=tau
	  ltau=ltau.new
	  tau=exp(ltau)
	  
	  if(max(abs(last.tau-tau)) / max(abs(last.tau)+abs(tau)) < control$nlminb$x.tol){
		break
	  }
	  n.nr=n.nr+1L
	  last.func=this.func
	  if(n.nr >= control$nlminb$iter.max){
		warning(paste(control$nlminb$iter.max, 'iterations reached'))
		break
	  }
	}    
  }else stop('Method not implemented')
  

  bd.idx=which(tau<boundary.eps)
  if(FALSE){
	  if(length(bd.idx)>0L){#browser()
		preprocPREML(tau)
		eval(preprocScore)
		all.score=score()
		bd.idx=bd.idx[ all.score[bd.idx] < -boundary.eps ]
		
		if(length(bd.idx)==length(tau)) {
		  tau=rep(0,length(tau))
		}else if(length(bd.idx)>0L) {
		  bd.idx0= bd.idx [ all.score[bd.idx]==0 ]
		  if(length(bd.idx0) != length(bd.idx)){
			new.fit=Recall(y, rep(0,length(y)), , k[-bd.idx], information, method, boundary.eps, conv, nStepHalving, control$nlminb$iter.max, starts=tau[-bd.idx], plot.it, verbose)
			if(new.fit$PREML >= PREML()){
			  tau[-bd.idx]=new.fit$parms
			  tau[bd.idx]=0
			}
		  }
		}
	  }
  }

  if((nearZero = length(bd.idx))>0L){  # check boundary
	bd.idx.all= if(nearZero == 1) list(matrix(bd.idx)) else lapply(seq_len(nearZero), combn, x=bd.idx)
	cur.obj=obj(tau)
	tau.bak = tau
	for(i.0 in bd.idx.all){
		if(nrow(i.0)<nK){
			for(j.0 in seq_len(ncol(i.0))){
				this.0 = i.0[,j.0]
				setZero=function(tau){tau0=tau.bak; tau0[this.0]=0; tau0[-this.0]=tau; tau0}  # this.0
				objNeg0=function(tau) objNeg(setZero(tau))
				gradNeg0=function(tau) gradNeg(setZero(tau))[-this.0]
				hessNeg0=function(tau) hessNeg(setZero(tau))[-this.0, -this.0, drop=FALSE]
				nlminb.fit=nlminb(tau[-this.0], objNeg0, gradNeg0, hessNeg0, lower=rep(0,nK-length(this.0)), control=nlminb.control)
				if( -nlminb.fit$objective > cur.obj) {tau = setZero(nlminb.fit$par); cur.obj=-nlminb.fit$objective}
			}
		}else{
			if( (tmp = obj(rep(0, nK))) > cur.obj){ tau = rep(0, nK); cur.obj = tmp}
		}
	}
  }
}
	
	if(length(candidates) > 0){
		curObj = -Inf; ans.idx=NA_integer_
		for(i in seq_along(candidates)){
			preprocPREML(candidates[i])
			newObj = PREML()
			if(newObj > curObj){
				curObj=newObj
				ans.idx = i
			}
		}
		tau = candidates[ans.idx]
	}else{
		tau = 0
	}
	preprocPREML(tau)
	eval(preprocScore)
	eval(preprocInfo)
	updateNegGrad()
	updateNegHess()
	sigma2=crossprod(LIy)/n
	
  if(isTRUE(plot.it)){
	taus=exp(seq(log(1e-9), log(2*tau), length=500))
	taus=sort(unique(c(taus, seq(0, 2*tau, length=500))))
	objs=sapply(taus, obj)
	# x11()
	plot(taus, objs, xlab='tau', ylab='objective', type='o',main='Profiled residual log likelihood')
	abline(v=tau)
  }
  
  nm=names(K)
  if(is.null(nm)) nm = if(nK>0L) paste('varComp', seq_along(K), sep='.') else character(0L)
  ans=list(
	## varComp.fit specific block
	parms=structure(tau, names=nm, candidates = candidates),
	gradients=structure(-negGrad, names=nm), 
	hessian=structure(hess(tau), dimnames=list(nm,nm)),
	sigma2=drop(sigma2), 
	varComps=structure(drop(tau*sigma2), names=nm),
	n.iter=n.nr, PREML=PREML(),
	X.Q2=Q2, 
	residual.contrast=y, 
	working.cor=structure(k, names=nm),
	
	## varComp common block
	# na.action=NULL,
	# offset = NULL,
	# contrasts=NULL,
	# xzlevels = NULL,
	# terms = NULL, 
	call=match.call(), 
	nobs = length(Y), 
	control=control, 
	random.labels = nm, 
	doFit= TRUE,
	# frame=if(isTRUE(keepXYK)) NULL else parent.frame(), 
	
	X=if(isTRUE(keepXYK)) X else NULL,
	# qrx = if(isTRUE(keepXYK)) qrx else NULL, 
	Y=if(isTRUE(keepXYK)) Y else NULL, 
	K=if(isTRUE(keepXYK)) K else NULL
  )
  class(ans)='varComp'
  ans
  
}
