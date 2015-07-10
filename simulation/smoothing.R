library(gmp); library(bigPolyRat); library(varComp)
library(splines); library(fda); library(nlme);

x0=-10:110 
x=0:100
n=length(x)
knots = seq(0,100,length=27)
int.knots=knots[seq(25)+1L]
nknots=length(knots)-2
X=cbind(1,x)
int.x = 1:99 = x

p=poly(x0, degree=18)[x0%in%x,]
    p = p*8 ## magic that makes the plots similar to those in the paper
    
f1=10*p[,1]-p[,7]+p[,12]
f2=10*p[,1]+2*p[,2]-p[,7]+5*p[,18]
f3=10*p[,1]+2*p[,2]+(100-x)^2*p[,18]/2000
#f4=10*p[,1]+30*p[,2]+5*p[,18]

set.seed(200)
err=rnorm(n)/2

y1=f1+err; y3=f2+err; y3=f3+err; y4=f2+sqrt(20)*err

op=par(mfrow=c(2,2)); plot(x,y1);plot(x,y3);plot(x,y3);plot(x,y4);par(op)

if(TRUE){
    B.cs = bs(x, knots=1:99, intercept=TRUE)
    Omega.co=bsplinepen(create.bspline.basis(rangeval = range(x), breaks=x, norder=4), Lfdobj = 2)
    source("simulation/formOmega.R")
    range(Omega.co-formOmega(0,100,1:99))
    eig.O.cs=eigen(Omega.co, TRUE)
    X.cs=B.cs%*%eig.O.cs$vec[,-seq(n)] ## this space does not include intercept. Why??
    Z.cs=B.cs%*%sweep(eig.O.cs$vec[,seq(n)], 2, sqrt(eig.O.cs$val[seq(n)]), '/')
    
}else{
    library(lmeSplines)  ## this has rank n-2 Z matrix rather than Z?
    cs = smspline.v(x)
    X.cs = cs$Xs
    Z.cs = cs$Zs
    K.cs = tcrossprod(Z.cs)
    Z.cs.normalized = cholRoot(normalizeTrace(K.cs))
    eig.K.cs=eigen(K.cs, TRUE)
}

if(TRUE){
	B.co=bsplineS(x, breaks=knots, norder=4, nderiv=0, returnMatrix=FALSE)
	Omega.co=bsplinepen(create.bspline.basis(rangeval = range(knots), breaks=knots, norder=4), Lfdobj = 2)
	eig.O.co=eigen(Omega.co, symmetric = TRUE)
	X.co = B.co%*%eig.O.co$vectors[,-seq(nknots+2)]
	Z.co = B.co%*%sweep(eig.O.co$vectors[,seq(nknots+2)], 2, sqrt(eig.O.co$values[seq(nknots+2)]),'/')
	K.co=tcrossprod(Z.co)
}else{
	X.co=X
	Z.co = approx.Z(smspline(knots), knots, x)
	K.co = tcrossprod(Z.co)
	Z.co.normalized = cholRoot(normalizeTrace(K.co))
}

(fit1.cs=varComp(y1~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit1.cs.noNorm=varComp(y1~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE), normalizeTrace=FALSE))
(fit1.cs=varComp(y1~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit1.co=varComp(y1~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit1.co=varComp(y1~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))

#(fit1.cs.po=varComp(y1~x, varcov=list(K.cs), control=varComp.control('po', plot.it=TRUE))) #decartes=2
(fit1.co.po=varComp(y1~x, varcov=list(K.co), control=list('po', plot.it=TRUE)))
coef(fit1.co.po,'var.ratio')
fit1.co.preml=Vectorize(logLik(fit1.co.po,FUN=T))
op=par(mfrow=c(1,2)); curve(fit1.co.preml(x), 0, .1);curve(fit1.co.preml(x), .1, 2*coef(fit1.co.po,'tau')); par(op)

(fit2.cs=varComp(y2~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit2.cs=varComp(y2~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit2.co=varComp(y2~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit2.co=varComp(y2~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
#(fit2.cs.po=varComp(y2~x, varcov=list(K.cs), control=varComp.control('po', plot.it=TRUE))) #decartes=3
(fit2.co.po=varComp(y2~x, varcov=list(K.co), control=list('po', plot.it=TRUE)))
coef(fit2.co.po,'var.ratio')
fit2.co.preml=Vectorize(logLik(fit2.co.po,FUN=T))
op=par(mfrow=c(1,2)); curve(fit2.co.preml(x), 0, 5);curve(fit2.co.preml(x), 5, 2*coef(fit2.co.po,'tau')); par(op)

(fit3.cs=varComp(y3~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit3.cs=varComp(y3~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit3.co=varComp(y3~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit3.co=varComp(y3~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
#(fit3.cs.po=varComp(y3~x, varcov=list(K.cs), control=list('po', plot.it=TRUE)))
(fit3.co.po=varComp(y3~x, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
#coef(fit3.cs.po,'var.ratio')
coef(fit3.co.po,'var.ratio')
fit3.co.preml=Vectorize(logLik(fit3.co.po,FUN=T))
op=par(mfrow=c(1,2)); curve(fit3.co.preml(x), 0, 10);curve(fit3.co.preml(x), 10, 2*coef(fit3.co.po,'tau')); par(op)

(fit4.cs=varComp(y4~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit4.cs=varComp(y4~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit4.co=varComp(y4~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit4.co=varComp(y4~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
#(fit4.cs.po=varComp(y4~x, varcov=list(K.cs), control=list('po',plot.it=TRUE)))
(fit4.co.po=varComp(y4~x, varcov=list(K.co), control=list('po',plot.it=TRUE)))
coef(fit4.co.po,'var.ratio')
fit4.co.preml=Vectorize(logLik(fit4.co.po,FUN=T))
op=par(mfrow=c(1,2)); curve(fit4.co.preml(x), 0, 5);curve(fit4.co.preml(x), 5, 2*coef(fit4.co.po,'tau')); par(op)

save.image('smoothing.RData')
