##=============================================================================

##' @title Change format of CATA datasets to perform the package functions
##'
##' @usage
##' change_cata_format2(Data, nprod, nattr, nsub, nsess)
##'
##' @description
##' CATATIS and CLUSCATA operate on data where the blocks of variables are merged horizontally.
##' If you have a vertical format, you can use this function to change the format.
##' The first column must contain the sessions, the second the subjects, the third the products and the others the attributes. If you don't have sessions, then the first column must contain the subjects and the second the products.
##' Unlike change_cata_format function, you can enter data with sessions and/or mixed data in terms of products/subjects. However, you have to set columns to indicate this beforehand.

##' @param Data data frame or matrix. Correspond to your data
##'
##' @param nprod integer. Number of products
##'
##' @param nattr integer. Number of attributes
##'
##' @param nsub integer. Number of subjects.
##'
##' @param nsess integer. Number of sessions
##'
##'
##' @return The arranged data for CATATIS and CLUSCATA function and the subjects names in the correct order.
##'
##'
##'
##' @keywords CATA
##'
##' @examples
##'
##' #Vertical format with sessions
##' data("PB_tuna")
##' chang=change_cata_format2(PB_tuna, nprod= 6, nattr= 27, nsub = 12, nsess= 3)
##' res.cat2=catatis(Data= chang$Datafinal, nblo = 12, NameBlocks =  chang$NameSub)
##'
##' #Vertical format without sessions
##' Data=PB_tuna[1:66,2:30]
##' chang2=change_cata_format2(Data, nprod= 6, nattr= 27, nsub = 11, nsess= 1)
##' res.cat3=catatis(Data= chang2$Datafinal, nblo = 11, NameBlocks =  chang2$NameSub)
##' res.clu3=cluscata(Data= chang2$Datafinal, nblo = 11, NameBlocks =  chang2$NameSub)
##'
##' @seealso   \code{\link{catatis}}, \code{\link{cluscata}}, \code{\link{change_cata_format}}
##'
##' @export


change_cata_format2=function(Data, nprod, nattr, nsub, nsess)
{
  Matjuge=list()
  Datafinal=NULL
  if (nsess==1)
  {
    Data=cbind(rep(1, nrow(Data)), Data)
  }
  Data=as.data.frame(Data)
  Data[,1]=as.factor(Data[,1])
  Data[,2]=as.factor(Data[,2])
  Data[,3]=as.factor(Data[,3])

  if(ncol(Data)!=nattr+3)
  {
    stop("The number of columns in the Data must be nattr+3 (or nattr+2 if you don't have sessions)")
  }

  sessions=levels(Data[,1])
  if(nlevels(Data[,1])!=nsess)
  {
    stop("The number of sessions is not equivalent to the number of levels in first column of Data")
  }
  subjects=levels(Data[,2])
  if(nlevels(Data[,2])!=nsub)
  {
    if (nsess==1)
    {
      stop("The number of subjects is not equivalent to the number of levels in first column of Data")
    }else{
      stop("The number of subjects is not equivalent to the number of levels in second column of Data")
    }
  }
  # products=levels(Data[,3])
  if(nlevels(Data[,3])!=nprod)
  {
    if (nsess==1)
    {
      stop("The number of products is not equivalent to the number of levels in second column of Data")
    }else{
      stop("The number of products is not equivalent to the number of levels in third column of Data")
    }
  }


  for (i in 1:nsub)
  {
    Matjuge[[i]]=matrix(0, nprod, nattr)
    Datajuge=Data[Data[,2]==subjects[i],]
    cptsess=0
    for (s in 1:nsess)
    {
      Datajugerep=Datajuge[Datajuge[,1]==sessions[s],]
      Datajugerep=Datajugerep[order(Datajugerep[,3]),]
      Datajugerep=as.data.frame((Datajugerep[,-c(1:3)]))
      if (nrow(Datajugerep)==nprod)
      {
        Matjuge[[i]]=Matjuge[[i]]+Datajugerep
        cptsess=cptsess+1
      }
    }
    Matjuge[[i]]=Matjuge[[i]]/cptsess
    if( i==1)
    {
      Datafinal=Matjuge[[i]]
    }else
    {
      Datafinal=cbind(Datafinal, Matjuge[[i]])
    }

  }
  rownames(Datafinal)=unique(Datajuge[,3][order(Datajuge[,3])])
  return(list(Datafinal=Datafinal, NameSub=subjects))
}
