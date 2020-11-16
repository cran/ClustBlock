.crit_cataXj_fast=function(Xj, index, Stot){
  #Xj list

  n=nrow(Xj[[1]])
  nblo=length(Xj)
  nvar=ncol(Xj[[1]])

  S=Stot[index, index]



  if (length(index)>1)
  {
    ressvd=svd(S)
    lambda=ressvd$d[1]
  }else{
    lambda=1
  }

  Q=nblo-lambda

  return(Q)

}
