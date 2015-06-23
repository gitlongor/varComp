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
  diag(mat[-1,])=1
  mat[,N-1L]=-co[-N]/co[N]
  rts=eigen(mat, symm=FALSE, only=TRUE)$value ## actually eigen converts mat to a dense matrix
  rts=rts[order(abs(Im(rts)), abs(Re(rts)))]
  if(isTRUE(only.real)) Re(rts[Im(rts)==0]) else rts
}