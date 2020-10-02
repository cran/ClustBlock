
##=============================================================================


##' @title Performs the STATIS method on different blocks of quantitative variables
##'
##' @usage
##' statis(Data,Blocks,NameBlocks=NULL,Graph_obj=TRUE, Graph_weights=TRUE, scale=FALSE)
##'
##' @description
##' STATIS method on quantitative blocks. SUpplementary outputs are also computed
##'
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param Graph_obj logical. Show the graphical representation od the objects? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item RV: the RV matrix: a matrix with the RV coefficient between blocks of variables
##'          \item compromise: a matrix which is the compromise of the blocks (akin to a weighted average)
##'          \item weights: the weights associated with the blocks to build the compromise
##'          \item lambda: the first eigenvalue of the RV matrix
##'          \item overall error : the error for the STATIS criterion
##'          \item error_by_conf: the error by configuration (STATIS criterion)
##'          \item rv_with_compromise: the RV coefficient of each block with the compromise
##'          \item homogeneity: homogeneity of the blocks (in percentage)
##'          \item coord: the coordinates of each object
##'          \item eigenvalues: the eigenvalues of the svd decomposition
##'          \item inertia: the percentage of total variance explained by each axis
##'          \item error_by_obj: the error by object (STATIS criterion)
##'          \item scalefactors: the scaling factors of each block
##'          \item proj_config: the projection of each object of each configuration on the axes: presentation by configuration
##'          \item proj_objects: the projection of each object of each configuration on the axes: presentation by object
##'          }
##'
##'
##'
##'
##' @references
##' \itemize{
##' \item Lavit, C., Escoufier, Y., Sabatier, R., Traissac, P. (1994). The act (statis method). Computational 462 Statistics & Data Analysis, 18 (1), 97-119.\\
##' \item Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods.Application to sensometrics. Food Quality and Preference, in Press.
##'}
##' @importFrom grDevices dev.new rainbow
##' @importFrom graphics abline arrows barplot par plot points text title
##' @importFrom stats sd
##'
##'
##' @keywords quantitative
##'
##' @examples
##'
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  st=statis(Data=smoo, Blocks=rep(2,24),NameBlocks = NameBlocks)
##'  summary(st)
##'  #with variables scaling
##'  st2=statis(Data=smoo, Blocks=rep(2,24),NameBlocks = NameBlocks, Graph_weights=FALSE, scale=TRUE)
##'
##' @seealso   \code{\link{plot.statis}}, \code{\link{clustatis}}
##'
##' @export

##=============================================================================


statis=function(Data,Blocks,NameBlocks=NULL,Graph_obj=TRUE, Graph_weights=TRUE, scale=FALSE){
  # Data size (n, all the var)
  # Blocks vector given the nb of var of each block

  #initialisation
  n=nrow(Data)
  nblo=length(Blocks)
  J=rep(1:nblo , times =  Blocks )

  #rownames, colnames, NameBlocks
  if (is.null(NameBlocks)) NameBlocks=paste("B",1:nblo,sep="-")
  if(is.null(rownames(Data))) rownames(Data)=paste0("X", 1:nrow(Data))
  if(is.null(colnames(Data))) colnames(Data)=paste0("Y",1:ncol(Data))

  #parapet for numerical Data
  for (i in 1: ncol(Data))
  {
    if (is.numeric(Data[,i])==FALSE)
    {
      stop(paste("Error: the data must be numeric (column",i,")"))
    }
  }

  #parapet for number of objects
  if(n<3)
  {
    stop("At least 3 objects are required")
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
        stop(paste("Error: Column", i, "is constant"))
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


  X=scale(Data,center=TRUE,scale=scale) #X contains the centered (and scaled if necessary) data tables
  Wj=array(0,dim=c(n,n,nblo));             # association matrices



  #scale factors
  for(j in 1:nblo) {
    Xj=as.matrix(X[,J==j])
    Wj[,,j]=Xj%*%t(Xj)
  }
  muk=NULL
  for (i in 1:nblo)
  {
    mi=Wj[,,i]
    muk[i]=sqrt(sum(diag(t(mi)%*%mi)))
  }
  mu=mean(muk)
  facteurech=mu/muk

  #parapet for constant configurations
  if (sum(muk==0)>0)
  {
    pb=which(muk==0)
    stop(paste("Error: Configuration", pb, "is constant"))
  }



  # globally standardization of each data matrix
  # Computation of association matrices
  for(j in 1:nblo) {
    Xj=as.matrix(X[,J==j])
    Wj[,,j]=Xj%*%t(Xj)
    Wj[,,j]=Wj[,,j]/sqrt(sum(diag(Wj[,,j]%*%Wj[,,j])))  # standardisation so that ||Wj||=1
  }
  # RV matrix:
  RV=matrix(0,nblo,nblo)
  diag(RV)=rep(1,nblo)
  if(nblo>1)
  {
    for (i in 1:(nblo-1)) {
      for (j in (i+1):nblo) {
        RV[i,j]=sum(diag(Wj[,,i]%*%Wj[,,j]))
        RV[j,i]=RV[i,j]
      } }
  }


  #first eigenvector and eigenvalue of RV matrix
  ressvd=svd(RV)
  u=ressvd$u[,1]
  u=u*sign(u[1])
  lambda=ressvd$d[1]

  # the compromise W:
  W=matrix(0,n,n)
  for (j in 1:nblo) { W=W+(u[j]*Wj[,,j]) }


  # error computation
  dw=rep(0,nblo)
  erreur=matrix(0,n,nblo)
  normW=sum(diag(W%*%t(W)))
  for (j in 1:nblo) {
    a=Wj[,,j]-(u[j]*W) #difference
    dw[j]=sum(diag(a%*%t(a))) #error by block
    erreur[,j]=diag(a%*%t(a))
  }
  Q=sum(dw) #statis criterion

  #error by object
  obj=rep(0,n)
  for (i in 1:n)
  {
    obj[i]=sum(erreur[i,])
  }

  #coordinates
  e=svd(W)
  C=e$u%*%sqrt(diag(abs(e$d)))


  #projection of each object of each block
  configs=array(0,c(n,n,nblo))
  if(nblo>1)
  {
    for (l in 1:nblo)
    {
      configs[,,l]=Wj[,,l]%*%(e$u%*%diag(c(sqrt(1/e$d[-n]),0)))*sqrt(lambda)
    }
  }


  #compute the presentation by object
  objects=array(0,c(nblo,n,n))
  if (nblo>1)
  {
    for (r in 1:nblo)
    {
      for (i in 1:n)
      {
        for (j in 1:n)
        {
          objects[ r, j, i] = configs[ i, j, r]
        }
      }
    }
  }


    #inertia of axes
    pouriner=round(e$d/sum(e$d)*100,2)
    pouriner=pouriner[-length(pouriner)]
    names(pouriner)=paste("Dim", 1:(nrow(Data)-1))

    #Graphical representation
    if (Graph_obj==TRUE)
    {
      dev.new()
      par(xpd=FALSE)
      plot(C[,1],C[,2],type="n",lwd=5,pch=16,xlab=paste("Dim 1 (",pouriner[1],"%)"), ylab=paste("Dim 2 (",pouriner[2],"%)"),xlim=c(min(C[,1])-0.2,max(C[,1])+0.2),ylim=c(min(C[,2])-0.2,max(C[,2])+0.2))
      text(C[,1],C[,2],rownames(Data),col=rainbow(nrow(Data)))
      abline(h=0,v=0)
      title("STATIS")

      #projection of each object of each block
      dev.new()
      plot(C[,1],C[,2],type="n",lwd=5,pch=16,xlab=paste("Dim 1 (",pouriner[1],"%)"), ylab=paste("Dim 2 (",pouriner[2],"%)"),xlim=c(min(C[,1])-0.25,max(C[,1])+0.25),ylim=c(min(C[,2])-0.25,max(C[,2])+0.25))
      text(C[,1],C[,2],rownames(Data),col=rainbow(n),font=2)
      for (l in 1: length(Blocks))
      {
        points(configs[,1,l],configs[,2,l],col=rainbow(n),pch=20)
        arrows(C[,1],C[,2],configs[,1,l],configs[,2,l],col=rainbow(n),lwd=0.5,lty = 2,angle=0)
      }

      abline(h=0,v=0)
      title("STATIS")

    }






  #RV with the compromise
  rv=NULL
  W_k=as.matrix(W)
  normW=sqrt(sum(diag(W_k%*%W_k)))
  for ( i in 1:nblo)
  {
    W_i=Wj[,,i]
    rv=c(rv,sum(diag(W_i%*%W_k))/(normW))
  }

  #homogeneity
  hom=lambda/nblo

  #remove the last axis
  configs=configs[,1:(n-1),]
  objects=objects[,1:(n-1),]




  #names
  rownames(RV)=colnames(RV)=names(u)=names(dw)=names(rv)=names(facteurech)=NameBlocks
  rownames(W)=colnames(W)=rownames(C)=names(obj)=rownames(Data)
  C=C[,-ncol(C)]
  rownames(C)=rownames(Data)
  colnames(C)=paste("Dim",1:(n-1))
  vp=e$d[-length(e$d)]
  colnames(C)=names(vp)=paste("Dim", 1:ncol(C))

  if(nblo>1)
  {
    dimnames(objects)[[1]]=dimnames(configs)[[3]]=NameBlocks
    dimnames(configs)[[2]]=dimnames(objects)[[2]]=paste("Dim", 1:ncol(C))
    dimnames(configs)[[1]]=dimnames(objects)[[3]]=rownames(Data)
  }

  if(Graph_weights==TRUE)
  {
    dev.new()
    barplot(u)
    title(paste("Weights"))
  }
  homogeneity=round(hom,3)*100

  res=list(RV=round(RV,2),compromise=round(W,4),weights=round(u,5),lambda=round(lambda,3),
           overall_error=round(Q,2),error_by_conf=round(dw,2),rv_with_compromise=round(rv,2),
           homogeneity=homogeneity,coord=round(C,5), eigenvalues=round(vp,3), inertia=pouriner, error_by_obj=round(obj,2),
           scalefactors=round(facteurech,2), proj_config=round(configs,5),
           proj_objects=round(objects,5), param=list(nblo=nblo, n=n))

  class(res)="statis"


  return(res)

}
