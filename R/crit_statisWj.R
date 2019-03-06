.crit_statisWj=function(Wj){
#Wj list



  n=nrow(Wj[[1]])
  nblo=length(Wj)


  # % RV matrix:
  RV=matrix(0,nblo,nblo)
  diag(RV)=rep(1,nblo)
  if(nblo>1)
  {
    for (i in 1:(nblo-1)) {
      for (j in (i+1):nblo) {
        RV[i,j]=sum(diag(Wj[[i]]%*%Wj[[j]]))
        RV[j,i]=RV[i,j]
      } }
  }



  ressvd=.firstEig(RV)
  u=ressvd$u
  u=u*sign(u[1])
  lambda=ressvd$lambda
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
