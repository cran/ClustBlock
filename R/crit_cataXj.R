.crit_cataXj=function(Xj, index, Stot){
  #Xj list

  n=nrow(Xj[[1]])
  nblo=length(Xj)
  nvar=ncol(Xj[[1]])

  S=Stot[index, index]

  if (length(index)>1)
  {
    ressvd=svd(S)
    u=ressvd$u[,1]
    u=u*sign(u[1])
    lambda=ressvd$d[1]
  }else{
    u=1
    lambda=1
  }

  # the compromise C:
  C=matrix(0,n,nvar)
  for (j in 1:nblo) { C=C+(u[j]*Xj[[j]]) }


  Q=nblo-lambda

  return(list(S=S,C=C,alpha=u,lambda=lambda,Q=Q))

}
