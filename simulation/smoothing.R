library(gmp); library(bigPolyRat); library(varComp)
library(splines); library(fda)

x0=-10:110 
x=0:100
n=length(x)
knots = c(-1, seq(0,100,length=25), 101)
nknots=length(knots)-2
X=cbind(1,x)
brks = c(-1, x, 101)

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

(fit1.cs=varComp(y1~0+X.cs, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit1.cs=varComp(y1~0+X.cs, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit1.co=varComp(y1~0+X.co, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit1.co=varComp(y1~0+X.co, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit1.co.po=varComp(y1~0+X.co, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit1.co.po,'var.ratio')

(fit2.cs=varComp(y2~0+X.cs, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit2.cs=varComp(y2~0+X.cs, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit2.co=varComp(y2~0+X.co, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit2.co=varComp(y2~0+X.co, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit2.co.po=varComp(y2~0+X.co, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit2.co.po,'var.ratio')

(fit3.cs=varComp(y3~0+X.cs, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit3.cs=varComp(y3~0+X.cs, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit3.co=varComp(y3~0+X.co, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit3.co=varComp(y3~0+X.co, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit3.cs.po=varComp(y3~0+X.cs, varcov=list(K.cs), control=list('po',start=1e4, plot.it=TRUE)))
(fit3.co.po=varComp(y3~0+X.co, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit3.cs.po,'var.ratio')
coef(fit3.co.po,'var.ratio')

(fit4.cs=varComp(y4~0+X.cs, varcov=list(K.cs), control=list(start=1e-3, plot.it=TRUE)))
(fit4.cs=varComp(y4~0+X.cs, varcov=list(K.cs), control=list(start=1e4, plot.it=TRUE)))
(fit4.co=varComp(y4~0+X.co, varcov=list(K.co), control=list(start=1e-3, plot.it=TRUE)))
(fit4.co=varComp(y4~0+X.co, varcov=list(K.co), control=list(start=1e4, plot.it=TRUE)))
(fit4.co.po=varComp(y4~0+X.co, varcov=list(K.co), control=list('po',start=1e4, plot.it=TRUE)))
coef(fit4.co.po,'var.ratio')
