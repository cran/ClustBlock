##=============================================================================


##' @title Test the consistency of the panel in a CATA experiment
##'
##' @description
##' Permutation test on the agreement between subjects in a CATA experiment
##'
##' @usage
##' consistency_cata_panel(Data,nblo, nperm=100, alpha=0.05)
##'
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param nblo  numerical. Number of blocks (subjects).
##'
##' @param nperm numerical. How many permutations are required? Default: 100
##'
##' @param alpha numerical between 0 and 1. What is the threshold? Default: 0.05
##'
##'
##'
##'@return a list with:
##'         \itemize{
##'          \item answer: the answer of the test
##'          \item pval: pvalue of the test
##'          \item dis: distance between the homogeneity and the median of the permutations
##'          }
##'
##'
##'
##' @keywords CATA
##'
##' @references
##' Llobell, F., Giacalone, D., Labenne, A.,  Qannari, E.M. (2019).	Assessment of the agreement and cluster analysis of the respondents in a CATA experiment.	Food Quality and Preference, 77, 184-190.
##' Bonnet, L., Ferney, T., Riedel, T., Qannari, E.M., Llobell, F. (September 14, 2022) .Using CATA for sensory profiling: assessment of the panel performance. Eurosense, Turku, Finland.
##'
##' @examples
##'\donttest{
##'  data(straw)
##' #with all subjects
##' consistency_cata_panel(Data=straw, nblo=114)
##'}
##'
##' @seealso   \code{\link{consistency_cata}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
##'
##' @export


##=============================================================================
consistency_cata_panel=function(Data,nblo, nperm=100, alpha=0.05)
{
  nprod=nrow(Data)
  nattr=ncol(Data)/nblo
  n=nprod
  Blocks=rep(nattr,nblo)
  J=rep(1:nblo , times =  Blocks )# indicates which block each variable belongs to
  X=Data

  #parapet for numerical Data
  for (i in 1: ncol(Data))
  {
    if (is.numeric(Data[,i])==FALSE)
    {
      stop(paste("The data must be numeric (column",i,")"))
    }
  }

  # #parapet for binary Data
  # if ((sum(Data==0)+sum(Data==1))!=(dim(Data)[1]*dim(Data)[2]))
  # {
  #   stop("only binary Data is accepted (0 or 1)")
  # }

  #no NA
  if(sum(is.na(Data))>0)
  {
    print("NA detected:")
    tabna=which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }


  Xj=array(0,dim=c(n,nattr,nblo))  # array with all subjects matrices
  muk=NULL
  for(j in 1:nblo)
  {
    Aj=as.matrix(X[,J==j])
    normXj=sqrt(sum(diag(tcrossprod(Aj, Aj))))
    muk[j]=normXj
    if(normXj==0)
    {
      stop(paste("error: the subject", j, "has only 0 or only 1"))  #parapet for constant configurations
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


  #first eigenvector and eigenvalue of S matrix
  ressvd=svd(S)
  lambda=ressvd$d[1]
  hom=lambda/nblo

  #permutations
  hom_perm=NULL
  lambda_perm=NULL
  for (k in 1:nperm)
  {

    for(j in 1:nblo)
    {
      Aj_perm=Aj=as.matrix(Xj[,,j])
      tirage=sample(1:nrow(X))
      for (i in 1:nrow(X))
      {
        Aj_perm[i,]=Aj[tirage[i],]
      }

      Xj[,,j]=Aj_perm
    }


    # S matrix:
    S=matrix(0,nblo,nblo)
    diag(S)=rep(1,nblo)
    for (i in 1:(nblo-1)) {
      for (j in (i+1):nblo) {
        S[i,j]=sum(diag(tcrossprod(Xj[,,i],Xj[,,j])))
        S[j,i]=S[i,j]
      } }



    #first eigenvector and eigenvalue of S matrix
    ressvd=svd(S)
    lambda_perm[k]=ressvd$d[1]
    hom_perm[k]=lambda_perm[k]/nblo
  }
  pval=sum(hom_perm>=hom)/nperm
  if (pval<alpha)
  {
    consist="the panel is consistent"
  }else{
    consist="the panel is not consistent"
  }

  dis=hom-median(hom_perm)

  return(list(answer= consist, pval=pval, dis=dis))
}
