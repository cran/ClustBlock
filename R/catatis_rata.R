##=============================================================================
##' @title Perform the CATATIS method on different blocks from a RATA experiment
##'
##' @usage
##' catatis_rata(Data,nblo,NameBlocks=NULL, NameVar=NULL, Graph=TRUE, Graph_weights=TRUE,
##'  Test_weights=FALSE, nperm=100)
##'
##' @description
##' CATATIS method for RATA data. Additional outputs are also computed. Non-binary data are accepted and weights can be tested.
##'
##'
##' @param Data data frame or matrix where the blocks of variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param nblo integer. Number of blocks (subjects).
##'
##' @param NameBlocks string vector. Name of each block (subject). Length must be equal to the number of blocks. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param NameVar string vector. Name of each variable (attribute, the same names for each subject). Length must be equal to the number of attributes. If NULL, the colnames of the first block are taken. Default: NULL
##'
##' @param Graph logical. Show the graphical representation? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##' @param Test_weights logical. Should the the weights be tested? Default: FALSE
##'
##' @param nperm integer. Number of permutation for the weight tests. Default: 100
##'
##' @return a list with:
##'         \itemize{
##'          \item S: the S matrix: a matrix with the similarity coefficient among the subjects
##'          \item compromise: a matrix which is the compromise of the subjects (akin to a weighted average)
##'          \item weights: the weights associated with the subjects to build the compromise
##'          \item weights_tests: the weights tests results
##'          \item lambda:  the first eigenvalue of the S matrix
##'          \item overall error: the error for the CATATIS criterion
##'          \item error_by_sub: the error by subject (CATATIS criterion)
##'          \item error_by_prod: the error by product (CATATIS criterion)
##'          \item s_with_compromise: the similarity coefficient of each subject with the compromise
##'          \item homogeneity: homogeneity of the subjects (in percentage)
##'          \item CA: the results of correspondence analysis performed on the compromise dataset
##'          \item eigenvalues: the eigenvalues associated to the correspondence analysis
##'          \item inertia: the percentage of total variance explained by each axis of the CA
##'          \item scalefactors: the scaling factors of each subject
##'          \item param: parameters called
##'          }
##'
##'
##'
##' @keywords RATA
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2019). A new approach for the analysis of data and the clustering of subjects in a CATA experiment. Food Quality and Preference, 72, 31-39.\cr
##' Bonnet, L., Ferney, T., Riedel, T., Qannari, E.M., Llobell, F. (September 14, 2022) .Using CATA for sensory profiling: assessment of the panel performance. Eurosense, Turku, Finland.\cr
##' Bonnet, L., Llobell, F., Qannari, E.M. (Pangborn 2023). Assessment of the panel performance in a RATA experiment.
##'
##'
##'
##' @examples
##' #RATA data with session
##' data(RATAchoc)
##' chang2=change_cata_format2(RATAchoc, nprod= 12, nattr= 13, nsub = 9, nsess= 3)
##' res.cat4=catatis_rata(Data= chang2$Datafinal, nblo = 9, NameBlocks =  chang2$NameSub)
##' summary(res.cat4)
##'
##' #RATA data without session
##' Data=RATAchoc[1:108,2:16]
##' chang2=change_cata_format2(Data, nprod= 12, nattr= 13, nsub = 9, nsess = 1)
##' res.cat5=catatis_rata(Data= chang2$Datafinal, nblo = 9, NameBlocks =  chang2$NameSub)
##' summary(res.cat5)
##' graphics.off()
##'
##' @seealso   \code{\link{catatis}}, \code{\link{plot.catatis}}, \code{\link{summary.catatis}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
##'
##' @export

##=============================================================================


catatis_rata=function(Data,nblo,NameBlocks=NULL, NameVar=NULL, Graph=TRUE, Graph_weights=TRUE, Test_weights=FALSE, nperm=100){

  #initialisation
  n=nrow(Data)
  p=ncol(Data)
  nvar=p/nblo
  #parapet for nblo
  if (as.integer(nvar)!=nvar)
  {
    stop("number of columns modulo nblo != 0")
  }
  Blocks=rep(nvar,nblo)
  J=rep(1:nblo , times =  Blocks )# indicates which block each variable belongs to


  #rownames, colnames, NameBlocks
  if (is.null(NameBlocks)) NameBlocks=paste("S",1:nblo,sep="-")
  if(is.null(rownames(Data))) rownames(Data)=paste0("X", 1:nrow(Data))
  if(is.null(colnames(Data))) colnames(Data)=rep(paste0("Y",1:nvar), nblo)

  X=Data


  #Parapet for NameBlocks
  if(length(NameBlocks)!=nblo)
  {
    stop("Error with the length of NameBlocks")
  }

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
    stop("At least 3 products are required")
  }

  #parapet for number of blocks
  if(nblo<2)
  {
    stop("At least 2 subjects are required")
  }

  #parapet for number of attributes
  if(nvar<3)
  {
    stop("At least 3 attributes are required")
  }

  #no NA
  if(sum(is.na(Data))>0)
  {
    print("NA detected:")
    tabna=which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }



  Xj=array(0,dim=c(n,nvar,nblo))  # array with all subjects matrices
  muk=NULL
  for(j in 1:nblo)
  {
    Aj=as.matrix(X[,J==j])
    normXj=sqrt(sum(diag(tcrossprod(Aj,Aj))))
    muk[j]=normXj
    if(normXj==0)
    {
      stop(paste("error: the subject",NameBlocks[j], "has only 0"))  #parapet for null configurations
    }
    Xj[,,j]=Aj/normXj #standardization
  }

  #scale factors
  mu=mean(muk)
  facteurech=mu/muk

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
  u=ressvd$u[,1]
  u=u*sign(u[1])
  lambda=ressvd$d[1]
  hom=lambda/sum(diag(S))

  # the compromise C:
  C=matrix(0,n,nvar)
  for (j in 1:nblo) { C=C+(u[j]*Xj[,,j]) }


  # the sum of distances between the weighted scalar product matrices and the consensus
  dw=rep(0,nblo)
  erreur=matrix(0,n,nblo)
  for (j in 1:nblo) {
    a=Xj[,,j]-(u[j]*C)
    dw[j]=sum(diag(tcrossprod(a)))
    erreur[,j]=diag(tcrossprod(a))
  }
  Q=sum(dw) #catatis criterion

  #error by object
  obj=rep(0,n)
  for (i in 1:n)
  {
    obj[i]=sum(erreur[i,])
  }

  #names of the compromise
  rownames(C)=names(obj)=rownames(Data)
  if (is.null(NameVar)==TRUE)
  {
    colnames(C)=colnames(Data)[1:nvar]
  }else{
    colnames(C)=NameVar
  }


  #s with compromise
  normC=sqrt(sum(diag(tcrossprod(C))))
  s=NULL
  for (i in 1:nblo)
  {
    s=c(s,sum(diag(tcrossprod(Xj[,,i],C)))/normC)
  }

  #CA
  compromis=C
  colomnnull=NULL
  for (l in 1:ncol(compromis))
  {
    if (sum(compromis[,l])==0)
    {
      colomnnull=c(colomnnull,l)
    }
  }
  rownull=NULL
  for (l in 1:nrow(compromis))
  {
    if (sum(compromis[l,])==0)
    {
      rownull=c(rownull,l)
    }
  }
  compromis2=compromis
  if(length(colomnnull)>0)
  {
    compromis2=compromis[,-colomnnull]
    warning("No block  has a 1 for the variable(s):  ", paste(colnames(compromis)[colomnnull], collapse=","))
  }
  if(length(rownull)>0)
  {
    compromis2=compromis2[-rownull,]
    warning("No block  has a 1 for the product(s):  ", paste(rownames(compromis)[rownull], collapse=","))
  }
  e=CA(compromis2,graph=FALSE)

  #inertia of axes
  pouriner=round(e$eig[,2],2)
  eigenvalues=round(e$eig[,1],4)

  #graphs
  if (Graph==TRUE)
  {
    dev.new()
    barplot(eigenvalues, col="blue", main="Eigenvalues")
    dev.new()
    options(ggrepel.max.overlaps = Inf)
    print(plot.CA(e,title=paste("CATATIS")))
  }

  #weight tests
  if (Test_weights==TRUE)
  {
    Weightsperm=NULL
    nperm=100
    tirage=sample(1:nblo, nperm, replace = TRUE)
    den=sum(diag(t(C)%*%C))
    for (i in 1:nperm)
    {
      Mat=Xj[,,tirage[i]]
      tirage2=sample(1:n)
      Matpermut=as.matrix(Mat[tirage2,])
      Weightsperm=c(Weightsperm, sum(diag(t(C)%*%Matpermut))/den)
    }
    valcrit=t.test(Weightsperm, alternative = "l")$conf.int[2]
    meanperm=mean(Weightsperm)

    Matrep=matrix(0, nrow=nblo, ncol=5)
    Matrep=as.data.frame(Matrep)
    Matrep[,1]=round(u,3)
    Matrep[,2]=rep(round(meanperm,3), nblo)
    Matrep[,3]=rep(round(valcrit,3), nblo)
    for (l in 1:nblo)
    {
      if (u[l]>valcrit)
      {
        Matrep[l,4]=round(u[l]-valcrit, 3)
        Matrep[l,5]="Yes"
      }else{
        Matrep[l,4]=NA
        Matrep[l,5]="No"
      }
    }
    colnames(Matrep)=c("Weights", "Random mean", "Critical value", "Dist with CV", "Signif")
    rownames(Matrep)= NameBlocks
  }else{
    Matrep=NULL
  }



  names(u)=names(s)=rownames(S)=colnames(S)= names(dw)=names(muk)=names(facteurech)=NameBlocks


  if(Graph_weights==TRUE)
  {
    dev.new()
    barplot(u)
    title(paste("Weights"))
  }
  homogeneity=round(hom,3)*100

  #results
  res=list(S=round(S,2),compromise=round(C,2),weights=round(u,5), weights_tests= Matrep, lambda=round(lambda,2),overall_error=round(Q,2),
           error_by_sub=round(dw,2), error_by_prod=round(obj,2), s_with_compromise=round(s,2), homogeneity=homogeneity, CA=e, eigenvalues=eigenvalues,
           inertia=pouriner, scalefactors=round(facteurech,2), param=list(n=n, nblo=nblo, nvar=nvar))
  class(res)="catatis"


  return(res)

}
