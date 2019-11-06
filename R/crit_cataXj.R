.crit_cataXj=function(Xj, index, Stot){
  #Xj list

  n=nrow(Xj[[1]])
  nblo=length(Xj)
  nvar=ncol(Xj[[1]])

  S=Stot[index, index]



  if (length(index)>1)
  {
    ressvd=.firstEig(S)
    u=ressvd$u
    u=u*sign(u[1])
    lambda=ressvd$lambda
  }else{
    u=1
    lambda=1
  }

  # the compromise C:
  C=matrix(0,n,nvar)
  for (j in 1:nblo) { C=C+(u[j]*Xj[[j]]) }


  # the sum of distances between the weighted scalar product matrices
  # and the consensus
  dw=rep(0,nblo)
  normW=sum(diag(C%*%t(C)))
  for (j in 1:nblo) {
    a=Xj[[j]]-(u[j]*C)
    dw[j]=sum(diag(a%*%t(a)))
  }
  Q=sum(dw)

  return(list(S=S,C=C,alpha=u,lambda=lambda,Q=Q))

}
