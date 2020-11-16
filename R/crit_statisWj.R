.crit_statisWj=function(Wj, index, RVtot){
  #Wj list



  n=nrow(Wj[[1]])
  nblo=length(Wj)

  RV=RVtot[index, index]

  if (length(index)>1)
  {
    ressvd=svd(RV)
    u=ressvd$u[,1]
    u=u*sign(u[1])
    lambda=ressvd$d[1]
  }else{
    u=1
    lambda=1
  }

  # the compromise W:
  W=matrix(0,n,n)
  for (j in 1:nblo) { W=W+(u[j]*Wj[[j]]) }

  Q=nblo - lambda

  return(list(RV=RV,W=W,u=u,lambda=lambda,Q=Q))

}
