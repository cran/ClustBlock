.crit_statisWj=function(Wj, index, RVtot){
  #Wj list



  n=nrow(Wj[[1]])
  nblo=length(Wj)

  RV=RVtot[index, index]

  if (length(index)>1)
  {
    ressvd=.firstEig(RV)
    u=ressvd$u
    u=u*sign(u[1])
    lambda=ressvd$lambda
  }else{
    u=1
    lambda=1
  }

  # the compromise W:
  W=matrix(0,n,n)
  for (j in 1:nblo) { W=W+(u[j]*Wj[[j]]) }


  # the sum of distances between the weighted scalar product matrices
  # and the consensus
  dw=rep(0,nblo)
  normW=sum(diag(W%*%t(W)))
  for (j in 1:nblo) {
    a=Wj[[j]]-(u[j]*W)
    dw[j]=sum(diag(a%*%t(a)))
  }
  Q=sum(dw)

  return(list(RV=RV,W=W,u=u,lambda=lambda,Q=Q))

}
