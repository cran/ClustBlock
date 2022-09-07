##=============================================================================


##' @title Preprocessing for Just About Right Data
##'
##' @description
##' For JAR data, this preprocessing is needed.
##'
##' @usage
##' preprocess_JAR(Data,  nprod, nsub, levelsJAR=3, beta=0.1)
##'
##'
##' @param Data data frame where the first column is the Assessors, the second is the products and all other columns the JAR attributes with numbers (1 to 3 or 1 to 5, see levelsJAR)
##'
##' @param nprod integer. Number of products.
##'
##' @param nsub integer. Number of subjects.
##'
##' @param levelsJAR integer. 3 or 5 levels. If 5, the data will be transformed in 3 levels.
##'
##' @param beta numerical. Parameter for agreement between JAR and other answers. Between 0 and 0.5.
##'
##' @return A list with:
##'         \itemize{
##'          \item Datafinal: the Data transformed
##'          \item NameSub: the name of each subject in the right order
##'          }
##'
##'
##' @keywords JAR
##'
##' @references
##' Llobell, F., Vigneau, E. & Qannari, E. M. (September 14, 2022). Multivariate data analysis and clustering of subjects in a Just about right task. Eurosense, Turku, Finland.
##'
##'
##' @examples
##' data(cheese)
##' prepro=preprocess_JAR(cheese, nprod=8, nsub=72, levelsJAR=5)
##'
##' @seealso   \code{\link{catatis_jar}}, \code{\link{cluscata_jar}}, \code{\link{cluscata_kmeans_jar}}
##'
##' @export


##=============================================================================

preprocess_JAR= function(Data,  nprod, nsub, levelsJAR=3, beta=0.1)
{
  NameSub=unique(Data[,1])
  NameProd=unique(Data[,2])
  JAR=Data
  #subject by subject
  for (i in 1:nsub)
  {
    JAR[((i-1)*nprod+1):(i*nprod),]=Data[Data[,1]==NameSub[i],]
  }
  #put data in same order of products
  for (i in 1:nprod)
  {
    JAR[seq(i,nsub*nprod+(i-1),nprod),]=Data[Data[,2]==NameProd[i],]
  }
  JAR=JAR[,-c(1:2)]

  #Use only three levels
  if (levelsJAR!=3 && levelsJAR!=5)
  {
    stop("levelsJAR parameter must be 3 or 5")
  }
  if (levelsJAR==5)
  {
    for (i in 1:nrow(JAR))
    {
      for (j in 1: ncol(JAR))
        if (JAR[i,j]==5 || JAR[i,j]==4)
        {
          JAR[i,j]=3
        }
      else if (JAR[i,j]==2)
      {
        JAR[i,j]=1
      }
      else if (JAR[i,j]==3)
      {
        JAR[i,j]=2
      }
    }
  }


  #rename the levels
  NameAttr=colnames(JAR)
  for (i in 1:nrow(JAR))
  {
    for (j in 1: ncol(JAR))
      if (JAR[i,j]==1)
      {
        JAR[i,j]=paste0(NameAttr[j], "_too_little")
      }
    else if (JAR[i,j]==2)
    {
      JAR[i,j]=paste0(NameAttr[j], "_jar")
    }
    else if (JAR[i,j]==3)
    {
      JAR[i,j]=paste0(NameAttr[j], "_too_much")
    }
  }

  #transform data in CATA data
  CATAData=NULL
  for (j in 1:ncol(JAR))
  {
    CATAData=cbind(CATAData, tab.disjonctif(JAR[,j]))
  }

  #use beta
  for (j in seq(1, ncol(JAR)-2, 3))
  {
    for (i in 1:nrow(CATAData))
    {
      if (CATAData[i, j]==1)
      {
        CATAData[i, c(j+1, j+2)] =beta/2
        CATAData[i,j]=1-beta
      }
      else if (CATAData[i, j+1]==1)
      {
        CATAData[i, j] =beta
        CATAData[i,j+1]=1-beta
      }
      else if (CATAData[i, j+2]==1)
      {
        CATAData[i, j] =beta
        CATAData[i,j+2]=1-beta
      }
    }
  }

  #Pass in horizontal format
  JAR2=change_cata_format(CATAData, nprod, ncol(CATAData), nsub = nsub, 1, NameProds = NameProd, NameAttr = colnames(CATAData))
  return(list(Datafinal=JAR2, NameSub= NameSub))
}
