
##' @title Compute the indices to evaluate the quality of the cluster partition in multi-block context
##'
##' @description
##' Compute the Il index to evaluate the agreement between each block and the global partition (in sensory: agreement between each subject and the global partition)
##'
##' Compute the Jl index to evaluate if each block has a partition (in sensory: if each subject made a partition of products)
##'
##' @usage
##' indicesClusters(Data, Blocks, cut, NameBlocks=NULL, center=TRUE, scale=FALSE)
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data.
##'
##' @param cut numerical vector. The partition of the cluster analysis.
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param center logical. Should the data variables be centered? Default: TRUE. Please set to FALSE for a CATA experiment
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##'
##' @return
##'         \itemize{
##'          \item Il: the Il indices
##'          \item jl: the jl indicess
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
##'
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
##' @seealso   \code{\link{clustRowsOnStatisAxes}}, , \code{\link{ClusMB}}
##'
##' @export


##=============================================================================

indicesClusters = function(Data, Blocks, cut, NameBlocks=NULL, center=TRUE, scale=FALSE)
{
  nblo=length(Blocks)
  n=nrow(Data)
  if (is.null(NameBlocks)) NameBlocks=paste("B",1:nblo,sep="-")
  if(is.null(rownames(Data))) rownames(Data)=paste0("X", 1:nrow(Data))
  if(is.null(colnames(Data))) colnames(Data)=paste0("Y",1:ncol(Data))
  nclust=length(unique(cut))

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

  J=rep(1:nblo , times =  Blocks )      # indicates which block each variable belongs to

  X=NULL
  Xj=list()
  resj=list()
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
    d=dist(Xi)
    truc=hclust(d, method="ward.D2")
    resj[[i]]=cutree(truc, nclust)
  }


  ####indices
  ind=NULL
  iner_intra=NULL
  iner_inter=NULL
  iner_tot=NULL
  ind_gen=NULL
  iner_intra_gen=NULL
  iner_inter_gen=NULL
  clustgen=cut
  for (l in 1:nblo)
  {
    clust=resj[[l]]
    inertie_intra=0
    inertie_inter=0
    inertie_intra_gen=0
    inertie_inter_gen=0
    cgglob=apply(Xj[[l]],2, mean)
    inertie_totale=sum((as.numeric(t(Xj[[l]]))-rep(cgglob,n))**2)
    for(k in 1:nclust)
    {
      if (sum(clust==k)>1)
      {
        cg=apply(Xj[[l]][clust==k,],2, mean)
      }
      else{
        cg= Xj[[l]][clust==k,]
      }

      if (sum(clustgen==k)>1)
      {
        cggen=apply(Xj[[l]][clustgen==k,],2, mean)
      }
      else{
        cggen= Xj[[l]][clustgen==k,]
      }

      inertie_intra=sum((as.numeric(t(Xj[[l]][clust==k,]))-rep(cg,sum(clust==k)))**2)+inertie_intra
      inertie_inter=sum(clust==k)*sum((cgglob-cg)**2)+inertie_inter

      inertie_intra_gen=sum((as.numeric(t(Xj[[l]][clustgen==k,]))-rep(cggen,sum(clustgen==k)))**2)+inertie_intra_gen
      inertie_inter_gen=sum(clustgen==k)*sum((cgglob-cggen)**2)+inertie_inter_gen
    }
    iner_intra=c(iner_intra, inertie_intra)
    iner_inter=c(iner_inter, inertie_inter)
    iner_tot=c(iner_tot, inertie_totale)
    ind=c(ind, inertie_inter/inertie_totale)

    iner_intra_gen=c(iner_intra_gen, inertie_intra_gen)
    iner_inter_gen=c(iner_inter_gen, inertie_inter_gen)
    ind_gen=c(ind_gen, inertie_inter_gen/inertie_totale)
  }
  names(ind)= names(iner_intra)=NameBlocks
  dev.new()
  barplot(ind, col="blue", ylab="Jl index", xlab="Blocks", main="Jl indices")

  names(ind_gen)= names(iner_intra_gen)=NameBlocks
  dev.new()
  barplot(ind_gen, col="blue", ylab="Il index", xlab="Blocks", main="Il indices")

  return(list(indicesJl= ind, indicesIl= ind_gen))
}
