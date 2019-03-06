.crit_cataXj=function(Xj){
  #Xj list

  n=nrow(Xj[[1]])
  nblo=length(Xj)
  nvar=ncol(Xj[[1]])


  # % S matrix:
  S=matrix(0,nblo,nblo)
  diag(S)=rep(1,nblo)
  if(nblo>1)
  {
    for (i in 1:(nblo-1))
    {
      for (j in (i+1):nblo)
      {
        S[i,j]=sum(diag(Xj[[i]]%*%t(Xj[[j]])))
        S[j,i]=S[i,j]
      }
    }
  }

  ressvd=.firstEig(S)
  u=ressvd$u
  u=u*sign(u[1])
  lambda=ressvd$lambda

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
