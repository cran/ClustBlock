
##' @title Perform a cluster analysis of rows in a Multi-block context with the ClusMB method
##'
##' @description
##' Clustering of rows (products in sensory analysis) in a Multi-block context.
##' The hierarchical clustering is followed by a partitioning algorithm (consolidation).
##'
##' @usage
##' ClusMB(Data, Blocks, NameBlocks=NULL, scale=FALSE, center=TRUE,
##' nclust=NULL, gpmax=6)
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data.
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param center logical. Should the data variables be centered? Default: TRUE. Please set to FALSE for a CATA experiment
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##' @param nclust numerical. Number of clusters to consider. If NULL, the Hartigan index advice is taken.
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default: min(6, number of blocks -2)
##'
##'
##' @return
##'         \itemize{
##'          \item group: the clustering partition after consolidation.
##'          \item nbgH: Advised number of clusters per Hartigan index
##'          \item nbgCH: Advised number of clusters per Calinski-Harabasz index
##'          \item cutree_k: the partition obtained by cutting the dendrogram in K clusters (before consolidation).
##'          \item dend: The ClusMB dendrogram
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##' @keywords quantitative CATA RATA
##'
##' @references
##' Llobell, F., Qannari, E.M. (June 10, 2022). Cluster analysis in a multi-bloc setting. SMTDA, Athens, Greece.\cr
##' Llobell, F., Giacalone, D., Qannari, E. M. (Pangborn 2021). Cluster Analysis of products in CATA experiments.\cr
##' Paper submitted
##'
##' @examples
##'
##' #####projective mapping####
##' library(ClustBlock)
##' data(smoo)
##' res1=ClusMB(smoo, rep(2,24))
##' summary(res1)
##' indicesClusters(smoo, rep(2,24), res1$group)
##'
##' ####CATA####
##' data(fish)
##' Data=fish[1:66,2:30]
##' chang2=change_cata_format2(Data, nprod= 6, nattr= 27, nsub = 11, nsess= 1)
##' res2=ClusMB(Data= chang2$Datafinal, Blocks= rep(27, 11), center=FALSE)
##' indicesClusters(Data= chang2$Datafinal, Blocks= rep(27, 11),cut = res2$group, center=FALSE)
##'
##' @seealso   \code{\link{indicesClusters}}, \code{\link{summary.clusRows}} , \code{\link{clustRowsOnStatisAxes}}
##'
##' @export


##=============================================================================


ClusMB= function(Data, Blocks, NameBlocks=NULL, scale=FALSE, center=TRUE,
                 nclust=NULL, gpmax=6)
{
  nblo=length(Blocks)
  n=nrow(Data)
  if (is.null(NameBlocks)) NameBlocks=paste("B",1:nblo,sep="-")
  if(is.null(rownames(Data))) rownames(Data)=paste0("X", 1:nrow(Data))
  if(is.null(colnames(Data))) colnames(Data)=paste0("Y",1:ncol(Data))

  #parapet for numerical Data
  for (i in 1: ncol(Data))
  {
    if (is.numeric(Data[,i])==FALSE)
    {
      stop(paste("The data must be numeric (column",i,")"))
    }
  }

  #parapet for number of objects
  if(n<3)
  {
    stop("At least 3 objects are required")
  }

  #parapet for number of blocks
  if(nblo<2)
  {
    stop("At least 2 blocks are required")
  }


  #parapet for Blocks
  if(sum(Blocks)!=ncol(Data))
  {
    stop("Error with Blocks")
  }

  #Parapet for NameBlocks
  if(length(NameBlocks)!=nblo)
  {
    stop("Error with the length of NameBlocks")
  }


  #parapet for scale: no constant variable
  if(scale==TRUE)
  {
    for (i in 1:ncol(Data))
    {
      if (sd(Data[,i])==0)
      {
        stop(paste("Column", i, "is constant"))
      }
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


  #centering and scaling if necessary
  Data=scale(Data,center=center,scale=scale)

  J=rep(1:nblo , times = Blocks) # indicates which block each variable belongs to

  X=NULL
  Xj=list()
  for (i in 1:nblo)
  {
    Xi=as.matrix(Data[,J==i])
    normXi=sqrt(sum(diag(Xi%*%t(Xi))))
    if(normXi==0)
    {
      stop(paste("error: the block",NameBlocks[i], "is constant"))  #parapet for constant configurations
    }
    Xi=Xi/normXi #standardization
    X=cbind(X, Xi)
    Xj[[i]]=Xi
  }

  d=dist(X)
  resh=hclust(d, method="ward.D2")
  res2=resh
  res2$height=(resh$height/sqrt(2))**2
  dev.new()
  plot(res2, hang=-1,  main="ClusMB Dendrogram", axes=TRUE,ylab="Height", sub="", xlab="" )

  #number of clusters advised by hartigan
  criter=sort(cumsum(res2$height),decreasing = TRUE)
  H=NULL
  for (k in 1:min(gpmax,n-2))
  {
    H[k]=(criter[k]/criter[k+1]-1)*(n-k-1)
  }
  nbgroup_hart=which.max(H[-(min(gpmax,n-2))]-H[-1])+1
  cat(paste("Recommended number of clusters =", nbgroup_hart),"\n")

  #number of clusters advised by Calinsi Harabasz
  bgss=NULL
  for (i in 1:gpmax)
  {
    bgss[i]=nblo-criter[i]
  }

  CH=NULL
  CH[1]=0
  for (k in 2:gpmax)
  {
    num=(bgss[k]/(k-1))
    den=(criter[k])/(n-k)
    CH[k]=num/den
  }
  nbgroup_ch=which.max(CH)

  #take the number of clusters suggested by Hartigan if not privided
  if (is.null(nclust))
  {
    nclust=nbgroup_hart
  }

  cutree_k=list()
  for (K in 1:gpmax)
  {
    coupe=cutree(res2,K)
    cutree_k[[K]]=coupe
  }
  names(cutree_k)=paste0("partition", 1:gpmax)

  clust=cutree(res2, nclust)
  #consolidation
  centers=matrix(0, nclust, ncol(X))
  for (i in 1:nclust)
    if (sum(clust==i)>1)
    {
      centers[i,]=apply(X[cutree(res2,nclust)==i,],2, mean)
    }else{
      centers[i,]=X[cutree(res2,nclust)==i,]
    }

  res=kmeans(X, centers, iter.max = 100)

  resClusMB=list(group = res$cluster, nbgH=nbgroup_hart,
                   nbgCH=nbgroup_ch, cutree_k= cutree_k, dend=res2,
                 param=list(nblo=nblo, gpmax=gpmax, n=n),
                 type="ClusMB")
  class(resClusMB)="clusRows"

  return(resClusMB)
}
