library(gmp); library(bigPolyRat); library(varComp)
library(splines); library(fda); library(nlme)

x0=-10:110 
x=0:100
n=length(x)
knots = seq(0,100,length=25)
nknots=length(knots)-2
X=cbind(1,x)
brks = x

p=poly(x0, degree=18)[x0%in%x,]
    p = p*8 ## magic that makes the plots similar to those in the paper
    
f1=10*p[,1]-p[,7]+p[,12]
f2=10*p[,1]+2*p[,2]-p[,7]+5*p[,18]
f3=10*p[,1]+2*p[,2]+(100-x)^2*p[,18]/2000

set.seed(23049)
err=rnorm(n)

y1=f1+err; y3=f2+err; y3=f3+err; y4=f2+sqrt(20)*err

op=par(mfrow=c(2,2)); plot(x,y1);plot(x,y3);plot(x,y3);plot(x,y4);par(op)


B.cs=bsplineS(x, breaks=brks, norder=4, nderiv=0, returnMatrix=FALSE)
Omega.cs=bsplinepen(create.bspline.basis(rangeval = range(brks), breaks=brks, norder=4), Lfdobj = 2)
eig.O.cs=eigen(Omega.cs, symmetric = TRUE)
X.cs = B.cs%*%eig.O.cs$vectors[,-seq(n+2)]
Z.cs = B.cs%*%eig.O.cs$vectors[,seq(n+2)]/sqrt(eig.O.cs$values[seq(n+2)])
K.cs=tcrossprod(Z.cs)

B.co=bsplineS(x, breaks=knots, norder=4, nderiv=0, returnMatrix=FALSE)
Omega.co=bsplinepen(create.bspline.basis(rangeval = range(knots), breaks=knots, norder=4), Lfdobj = 2)
eig.O.co=eigen(Omega.co, symmetric = TRUE)
X.co = B.co%*%eig.O.co$vectors[,-seq(nknots+2)]
Z.co = B.co%*%eig.O.co$vectors[,seq(nknots+2)]/sqrt(eig.O.co$values[seq(nknots+2)])
K.co=tcrossprod(Z.co)

(fit1.cs=varComp(y1~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit1.cs.noNorm=varComp(y1~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE), normalizeTrace=FALSE))
(fit1.cs=varComp(y1~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit1.co=varComp(y1~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit1.co=varComp(y1~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit1.co.po=varComp(y1~x, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit1.co.po,'var.ratio')


spars=seq(0.5,.51,length=1000)
lambs=numeric(length(spars))
for(i in seq_along(spars))lambs[i]=smooth.spline(x, y1, all.knots=TRUE, spar=spars[i])$lambda
i.opt=which.min(abs(lambs-coef(f0,'var.ratio')))
plot(lambs); abline(h=coef(fit1.cs.noNorm,'var.ratio'))
cs.fit1=smooth.spline(x, y1, all.knots=TRUE, spar=spars[i.opt])

ones=rep(1,n); Z.cs.normalized=cholRoot(normalizeTrace(K.cs))
fit1.cs.lme=lme(y1~x, random=list(ones=pdIdent(~-1+Z.cs.normalized)))
VarCorr(fit1.cs.lme) ## identical fit
plot(x, y1); lines(x, fitted(fit1.cs.lme)); lines(x, fitted(cs.fit1), col=4) # Different:(

vcov.Y1.cs=vcov(fit1.cs, what='Y')
xb.cs1=model.matrix(fit1.cs,'X')%*%coef(fit1.cs)
ranef.cs1=coef(fit1.cs, 'varComp')[1]*crossprod(model.matrix(fit1.cs,'Z')[[1]], solve(vcov.Y1.cs, y1-xb.cs1))
Y1.blup=xb.cs1 + model.matrix(fit1.cs,'Z')[[1]]%*%ranef.cs1
lines(x, Y1.blup, col=2) # same as lme fit

B.coef = solve(crossprod(B.cs)+cs.fit1$lambda*Omega.cs, crossprod(B.cs,y1))
Y1.cs = B.cs%*%B.coef
lines(x, Y1.cs, col=3)
    
(fit2.cs=varComp(y2~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit2.cs=varComp(y2~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit2.co=varComp(y2~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit2.co=varComp(y2~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit2.co.po=varComp(y2~x, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit2.co.po,'var.ratio')

(fit3.cs=varComp(y3~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit3.cs=varComp(y3~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit3.co=varComp(y3~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit3.co=varComp(y3~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit3.cs.po=varComp(y3~x, varcov=list(K.cs), control=list('po',start=1e4, plot.it=TRUE)))
(fit3.co.po=varComp(y3~x, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit3.cs.po,'var.ratio')
coef(fit3.co.po,'var.ratio')

(fit4.cs=varComp(y4~x, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit4.cs=varComp(y4~x, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit4.co=varComp(y4~x, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit4.co=varComp(y4~x, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit4.co.po=varComp(y4~x, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit4.co.po,'var.ratio')
