.crit_cataMean=function(Xj){
  #Xj list

  n=nrow(Xj[[1]])
  nblo=length(Xj)
  nvar=ncol(Xj[[1]])

  # the compromise C:
  C=matrix(0,n,nvar)
  for (j in 1:nblo) {C=C+Xj[[j]]}
  C=C/nblo

  # the sum of distances between the weighted scalar product matrices
  # and the consensus
  dw=rep(0,nblo)
  for (j in 1:nblo) {
    a=Xj[[j]]-C
    dw[j]=sum(diag(a%*%t(a)))
  }
  Q=sum(dw)

  return(list(C=C,Q=Q))

}
