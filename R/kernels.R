
normalizeTrace=function(x){x=as.matrix(x); x/mean(diag(x), na.rm=TRUE)}

Minkowski=function(x, p=1) 1-as.matrix(dist(x, method='minkowski', p=p)) * .5 / max(1, ncol(x))^(1/p) 

IBS=function(x)  1 - as.matrix(dist(x, method='manhattan') * .5 /max(1, ncol(x)) )  ## dist does scaling in the presence of missing values

Lin0=function(x) normalizeTrace(tcrossprod(x)/max(1,ncol(x)))
Quad1=function(x) normalizeTrace((base::tcrossprod(x)+1)^2/max(1,ncol(x)))
# Intxn2=function(K1, K2) normalizeTrace(K1*K2)
Polyk=function(x,c=0,d=1) normalizeTrace((base::tcrossprod(x)+c)^d)

cholRoot=function(x)
{
	R=suppressWarnings(chol(x, pivot=TRUE))
	oo <- order(attr(R, 'pivot'))
	rk = attr(R, 'rank')
	t(R[seq_len(rk), oo, drop=FALSE])
}


if(FALSE){
AM=function(x) #SPA3G::KERNEL(x, rep(1,max(1,ncol(x)))) ## this is not IBS kernel! The difference is 1 vs 1 comparison: IBS treat this as 2 (out of 2) but AM treat this as 2 (out of 4).
#below is copied from v1.0 of SPA3G (ARCHIVED on CRAN)
function (G, weight) 
{	weight=rep(1,max(1,ncol(x)))
    if (length(dim(G)) == 0) {
        size <- length(G)
        k <- matrix(1, size, size)
        for (i in 1:(size - 1)) {
            j <- seq(1, i, 1)
            remain <- G[-j]
            Ones <- matrix(1, length(remain), 1)
            leading <- Ones * G[i]
            D <- abs(remain - leading)
            AM <- D
            AM[AM == 0] <- 4
            AM[AM == 2] <- 0
            AM[AM == 1] <- 2
            AM[remain == 1 & leading == 1] <- 2
            k[i, (i + 1):size] <- k[(i + 1):size, i] <- AM * 
                weight/sum(4 * weight)
        }
    }
    if (length(dim(G)) > 0) {
        size <- nrow(G)
        k <- matrix(1, size, size)
        for (i in 1:(size - 1)) {
            j <- seq(1, i, 1)
            if (i < (size - 1)) {
                remain = as.matrix(G[-j, ])
            }
            if (i == (size - 1)) {
                remain <- t(as.matrix(G[-j, ]))
            }
            Ones <- matrix(1, nrow(remain), 1)
            leading <- Ones %*% G[i, ]
            D <- abs(remain - leading)
            AM <- as.matrix(D)
            AM[AM == 0] <- 4
            AM[AM == 2] <- 0
            AM[AM == 1] <- 2
            AM[remain == 1 & leading == 1] <- 2
            k[i, (i + 1):size] <- k[(i + 1):size, i] <- AM %*% 
                weight/sum(4 * weight)
        }
    }
    return(k)
}
am = function(x) cholRoot(AM(x))

}


ibs=function(x)	cholRoot(IBS(x))
lin0=function(x)	cholRoot(Lin0(x))
quad1=function(x) cholRoot(Quad1(x))
polyk=function(x,c=0,d=1) cholRoot(Polyk(x,c,d))
minkowski=function(x,p=1) cholRoot(Minkowski(x,p))

# fbr2=function(fixed,rdm) {
	# fixedX=model.matrix(~fixed)[,-1L, drop=FALSE]
	# ans = lapply(seq_len(ncol(fixedX)), function(i)fixedX[,i]*rdm)
	# ans=do.call('cbind', ans)
	# names = outer(colnames(fixedX), colnames(rdm), paste0, sep=':')
	# grp = paste0('.__fbr2[',seq_len(ncol(fixedX)), ']')[col(names)]
	# grpNames = paste0(grp, names, sep='_')
	# colnames(ans)=grpNames
	# ans
# }

#TODO: add more kernels, e.g., Gaussian
