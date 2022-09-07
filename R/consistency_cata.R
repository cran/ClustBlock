##=============================================================================


##' @title Test the consistency of each attribute in a CATA experiment
##'
##' @description
##' Permutation test on the agreement between subjects for each attribute in a CATA experiment
##'
##' @usage
##' consistency_cata(Data,nblo, nperm=100, alpha=0.05, printAttrTest=FALSE)
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
##' @param printAttrTest logical. Print the number of remaining attributes to be tested? Default: FALSE
##'
##'
##'
##'@return a list with:
##'         \itemize{
##'          \item consist: the consistent attributes
##'          \item no_consist: the inconsistent attributes
##'          \item pval: pvalue for each test
##'          }
##'
##'
##'
##' @keywords CATA
##'
##' @references
##' Llobell, F., Giacalone, D., Labenne, A.,  Qannari, E.M. (2019).	Assessment of the agreement and cluster analysis of the respondents in a CATA experiment.	Food Quality and Preference, 77, 184-190.
##'
##'
##' @examples
##'\donttest{
##'  data(straw)
##' #with only 40 subjects
##' consistency_cata(Data=straw[,1:(16*40)], nblo=40)
##' #with all subjects
##' consistency_cata(Data=straw, nblo=114, printAttrTest=TRUE)
##'}
##'
##' @seealso   \code{\link{consistency_cata_panel}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
##'
##' @export


##=============================================================================

consistency_cata=function(Data,nblo, nperm=100, alpha=0.05, printAttrTest=FALSE)
{
  nprod=nrow(Data)
  nattr=ncol(Data)/nblo
  n=nprod
  Blocks=rep(nattr,nblo)
  J=rep(1:nblo , times =  Blocks )# indicates which block each variable belongs to
  if(is.null(colnames(Data))) colnames(Data)=rep(paste0("A",1:Blocks[1]), nblo)

  #parapet for numerical Data
  for (i in 1: ncol(Data))
  {
    if (is.numeric(Data[,i])==FALSE)
    {
      stop(paste("The data must be numeric (column",i,")"))
    }
  }

  #no NA
  if(sum(is.na(Data))>0)
  {
    print("NA detected:")
    tabna=which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }


  #prepare data
  Xi=array(0,dim=c(n,Blocks[1],nblo))  # array with all subjects matrices
  for(j in 1:nblo)
  {
    Xi[,,j]=as.matrix(Data[,J==j])
  }

  Tabi=array(0, dim=c(nprod, nblo, nattr))
  for (j in 1:nattr)
  {
    for (k in 1:nblo)
    {
      normk=sqrt(sum(diag(tcrossprod(Xi[,j,k],Xi[,j,k]))))
      if(normk>0)
      {
        Tabi[,k,j]=Xi[,j,k]/normk
      }
      else{
        Tabi[,k,j]=Xi[,j,k]
      }
    }
  }

  # S matrices:
  Sarray=array(0, dim=c(nblo,nblo,nattr))
  for (k in 1:nattr)
  {
    S=matrix(0,nblo,nblo)
    diag(S)=rep(1,nblo)
    for (i in 1:(nblo-1))
    {
      for (j in (i+1):nblo)
      {
        S[i,j]=sum(diag(tcrossprod(Tabi[,i,k],Tabi[,j,k])))
        S[j,i]=S[i,j]
      }
      Sarray[,,k]=S
    }
  }

  l1=NULL
  u1=NULL
  for (i in 1:nattr)
  {
    ei=eigen(Sarray[,,i])
    lambda1=ei$values[1]
    l1=c(l1, lambda1)
    u1=cbind(u1, abs(ei$vectors[,1]))
  }
  names(l1)=colnames(Data[,1:nattr])
  hom=(l1/nblo)*100

  # Permutations:
  res_perm2=matrix(0, nperm, nattr)
  for (k in 1:nattr)
  {
    lambdaun_perm=NULL
    eig=NULL
    Tabi_per=Tabi
    for (perm in 1:nperm)
    {
      for(j in 1:nblo)
      {
        Tabi_per[,j,k]=Tabi[sample(1:n),j,k]
      }

      S_perm=matrix(0,nblo,nblo)
      diag(S_perm)=rep(1,nblo)
      for (i in 1:(nblo-1))
      {
        for (j in (i+1):nblo)
        {
          S_perm[i,j]=sum(diag(tcrossprod(Tabi_per[,i,k],Tabi_per[,j,k])))
          S_perm[j,i]=S_perm[i,j]
        }
      }
      eig=c(eig,eigen(S_perm)$values[1])
    }
    res_perm2[,k]=eig
    if(printAttrTest==TRUE)
    {
      print(nattr-k)
    }
  }

  #results
  quant=NULL
  consist=NULL
  no_consist=NULL
  pval=NULL
  for (k in 1:nattr)
  {
    pval[k]=sum(l1[k]<res_perm2[,k])/nperm
    if(pval[k]>alpha)
    {
      no_consist=c(no_consist,colnames(Data)[k])
    }else{
      consist=c(consist,colnames(Data)[k])
    }
  }

  names(pval)=names(l1)

  return(list(consist=consist, no_consist=no_consist, pval=pval))
}
