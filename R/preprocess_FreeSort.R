##=============================================================================


##' @title Preprocessing for Free Sorting Data
##'
##' @description
##' For Free Sorting Data, this preprocessing is needed.
##'
##' @usage
##' preprocess_FreeSort(Data, NameSub=NULL)
##'
##'
##' @param Data data frame or matrix. Corresponds to all variables that contain subjects results. Each column corresponds to a subject and gives the groups to which the products (rows) are assigned
##'
##' @param NameSub string vector. Name of each subject. Length must be equal to the number of clumn of the Data. If NULL, the names are S1,...Sm. Default: NULL
##'
##'
##' @return A list with:
##'         \itemize{
##'          \item new_Data: the Data transformed
##'          \item Blocks: the number of groups for each subject
##'          \item NameBlocks: the name of each subject
##'          }
##'
##'
##' @keywords FreeSorting
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods. Application to sensometrics. Food Quality and Preference, in Press.
##'
##'
##' @importFrom  FactoMineR tab.disjonctif
##'
##'
##' @examples
##'data(choc)
##'prepro=preprocess_FreeSort(choc)
##'
##' @seealso   \code{\link{clustatis}}, \code{\link{clustatis_FreeSort}}
##'
##' @export


##=============================================================================


preprocess_FreeSort=function(Data, NameSub=NULL)
{
  p=ncol(Data)
  n=nrow(Data)

  if(is.null(colnames(Data))) colnames(Data)=paste0("S",1:ncol(Data))
  if (is.null(NameSub)) NameSub=paste("S",1:p)

  #parapet for Namesub
  if(length(NameSub)!=p)
  {
    stop("Namesub must have the length of ncol(Data)")
  }


  Yi=list()
  for (i in 1:p)
  {
    Data[,i]=factor(Data[,i])
    Yi[[i]]=as.matrix(tab.disjonctif(Data[,i]))
    a=apply(Yi[[i]],2,mean)
    Yi[[i]]=Yi[[i]]%*%diag(1/sqrt(a))#standardization
    colnames(Yi[[i]])=paste0(NameSub[i],"g",1:ncol(Yi[[i]]))
  }

  dat=NULL
  for (i in 1:p)
  {
    dat=cbind(dat,Yi[[i]])
  }
  rownames(dat)=rownames(Data)
  Blocks=rep(0,p)
  for (i in 1:p)
  {
    Blocks[i]=nlevels(as.factor(Data[,i]))
  }


  return(list(new_Data=dat, Blocks=Blocks,NameBlocks=NameSub))
}
