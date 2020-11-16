.vp_s=function(Data,nblo, scale=FALSE){
  # Data size (n, all the var)
  # Blocks vector given the nb of var of each block

  #initialisation
  n=nrow(Data)
  p=ncol(Data)
  nvar=p/nblo
  Blocks=rep(nvar,nblo)
  J=rep(1:nblo , times =  Blocks )# indicates which block each variable belongs to
  X=Data

  Xj=array(0,dim=c(n,Blocks[1],nblo))  # array with all subjects matrices
  muk=NULL
  for(j in 1:nblo)
  {
    Aj=as.matrix(X[,J==j])
    normXj=sqrt(sum(Aj==1))
    muk[j]=normXj
    if(normXj==0)
    {
      stop(paste("error: the subject",j, "has only 0"))  #parapet for constant configurations
    }
    Xj[,,j]=Aj/normXj #standardization
  }

  # S matrix:
  S=matrix(0,nblo,nblo)
  diag(S)=rep(1,nblo)
  for (i in 1:(nblo-1)) {
    for (j in (i+1):nblo) {
      S[i,j]=sum(diag(tcrossprod(Xj[,,i],Xj[,,j])))
      S[j,i]=S[i,j]
    } }


  #eigenvalues of S matrix
  ressvd=svd(S)
  lambda=ressvd$d
  return(lambda[-nblo])
}
