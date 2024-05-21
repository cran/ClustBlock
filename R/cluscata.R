##=============================================================================


##' @title Perform a cluster analysis of subjects from a CATA experiment
##'
##' @description
##' Clustering of subjects (blocks) from a CATA experiment. Each cluster of blocks is associated with a compromise
##' computed by the CATATIS method. The hierarchical clustering is followed by a partitioning algorithm (consolidation).
##' Non-binary data are accepted.
##'
##' @usage
##'cluscata(Data, nblo, NameBlocks=NULL, NameVar=NULL, Noise_cluster=FALSE,
##'         Itermax=30, Graph_dend=TRUE, Graph_bar=TRUE, printlevel=FALSE,
##'         gpmax=min(6, nblo-2), rhoparam=NULL, Testonlyoneclust=FALSE, alpha=0.05,
##'         nperm=50, Warnings=FALSE)
##'
##' @param Data data frame or matrix where the blocks of binary variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param nblo  numerical. Number of blocks (subjects).
##'
##' @param NameBlocks string vector. Name of each block (subject). Length must be equal to the number of blocks. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param NameVar string vector. Name of each variable (attribute, the same names for each subject). Length must be equal to the number of attributes. If NULL, the colnames of the first block are taken. Default: NULL
##'
##' @param Noise_cluster logical. Should a noise cluster be computed? Default: FALSE
##'
##' @param Itermax numerical. Maximum of iteration for the partitioning algorithm. Default:30
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Default: TRUE
##'
##' @param printlevel logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default: min(6, nblo-2)
##'
##' @param rhoparam numerical. What is the threshold for the noise cluster? Between 0 and 1, high value can imply lot of blocks set aside. If NULL, automatic threshold is computed.
##'
##' @param Testonlyoneclust logical. Test if there is more than one cluster? Default: FALSE
##'
##' @param alpha numerical between 0 and 1. What is the threshold to test if there is more than one cluster? Default: 0.05
##'
##' @param nperm numerical. How many permutations are required to test if there is more than one cluster? Default: 50
##'
##' @param Warnings logical. Display warnings about the fact that none of the subjects in some clusters checked an attribute or product? Default: FALSE
##'
##'
##' @return Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
##'         \itemize{
##'          \item group: the clustering partition after consolidation. If Noise_cluster=TRUE, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold for the noise cluster
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item s_with_compromise: similarity coefficient of each subject with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item compromise: the compromise of each cluster
##'          \item CA: list. the correspondance analysis results on each cluster compromise (coordinates, contributions...)
##'          \item inertia: percentage of total variance explained by each axis of the CA for each cluster
##'          \item s_all_cluster: the similarity coefficient between each subject and each cluster compromise
##'          \item criterion: the CLUSCATA criterion error
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'          There is also at the end of the list:
##'          \itemize{
##'          \item dend: The CLUSCATA dendrogram
##'          \item cutree_k: the partition obtained by cutting the dendrogram in K clusters (before consolidation).
##'          \item overall_homogeneity_ng: percentage of overall homogeneity by number of clusters before consolidation (and after if there is no noise cluster)
##'          \item diff_crit_ng: variation of criterion when a merging is done before consolidation (and after if there is no noise cluster)
##'          \item test_one_cluster: decision and pvalue to know if there is more than one cluster
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##' @keywords CATA
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2019). A new approach for the analysis of data and the clustering of subjects in a CATA experiment. Food Quality and Preference, 72, 31-39.\cr
##' Llobell, F., Giacalone, D., Labenne, A.,  Qannari, E.M. (2019).	Assessment of the agreement and cluster analysis of the respondents in a CATA experiment.	Food Quality and Preference, 77, 184-190.
##'
##' @examples
##' \donttest{
##' data(straw)
##' #with 40 subjects
##' res=cluscata(Data=straw[,1:(16*40)], nblo=40)
##' #plot(res, ngroups=3, Graph_dend=FALSE)
##' summary(res, ngroups=3)
##' #With noise cluster
##' res2=cluscata(Data=straw[,1:(16*40)], nblo=40, Noise_cluster=TRUE,
##' Graph_dend=FALSE, Graph_bar=FALSE)
##' #With noise cluster and defined rho threshold
##' #(high threshold for this example, you can put low threshold
##' #(ex: 0.2 or 0.3) to avoid set aside lot of respondents)
##' res3=cluscata(Data=straw[,1:(16*40)], nblo=40, Noise_cluster=TRUE,
##' Graph_dend=FALSE, Graph_bar=FALSE, rhoparam=0.6)
##' #with all subjects
##' res=cluscata(Data=straw, nblo=114, printlevel=TRUE)
##'
##'
##' #Vertical format
##' data("fish")
##' Data=fish[1:66,2:30]
##' chang2=change_cata_format2(Data, nprod= 6, nattr= 27, nsub = 11, nsess= 1)
##' res3=cluscata(Data= chang2$Datafinal, nblo = 11, NameBlocks =  chang2$NameSub)
##' }
##'
##' @seealso   \code{\link{plot.cluscata}}, \code{\link{summary.cluscata}} , \code{\link{catatis}}, \code{\link{cluscata_kmeans}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
##'
##' @export


##=============================================================================




cluscata=function(Data, nblo, NameBlocks=NULL, NameVar=NULL, Noise_cluster=FALSE, Itermax=30,
                  Graph_dend=TRUE, Graph_bar=TRUE, printlevel=FALSE,
                  gpmax=min(6, nblo-2), rhoparam=NULL,
                  Testonlyoneclust=FALSE, alpha=0.05, nperm=50, Warnings=FALSE){

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
  if(nblo<4)
  {
    stop("At least 4 subjects are required")
  }

  #parapet for number of attributes
  if(nvar<3)
  {
    stop("At least 3 attributes are required")
  }



  #parapet for gpmax
  if (gpmax>(nblo-1))
  {
    stop(paste("gpmax > number of blocks-1"))
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
  Xi=array(0,dim=c(n,nvar,nblo))  # array with all subjects matrices
  muk=NULL
  for(j in 1:nblo)
  {
    Aj=as.matrix(X[,J==j])
    normXi=sqrt(sum(diag(tcrossprod(Aj,Aj))))
    muk[j]=normXi
    if(normXi==0)
    {
      stop(paste("error: the subject",NameBlocks[j], "has only 0"))
    }
    Xi[,,j]=Aj/normXi #standardization
  }


  #only one cluster?
  if (Testonlyoneclust==TRUE)
  {
    testonecluster=.one_cluster_or_more_cata(Data, nblo, nperm = nperm, alpha=alpha)
    if(testonecluster$decision==TRUE)
    {
      testonecluster$decision="Only one cluster can be considered"
    }else{
      testonecluster$decision="Clustering is necessary"
    }
  }else{
    testonecluster=list(decision="Untested", pvalue="Untested")
  }

  # S matrix:
  S=matrix(0,nblo,nblo)
  diag(S)=rep(1,nblo)
  for (i in 1:(nblo-1)) {
    for (j in (i+1):nblo) {
      S[i,j]=sum(diag(tcrossprod(Xi[,,i],Xi[,,j])))
      S[j,i]=S[i,j]
    } }

  # criterion when each data table forms a group itself: 0
  crit=rep(0,nblo)
  Q_current=0
  critS=rep(1,nblo)
  S_current=nblo

  cc=1:nblo
  # the names of the clusters are 1 to 'number of clusters'
  code=1:nblo;
  # the names of the clusters are taken from the HCA
  cvar=1:nblo
  # when a new cluster is formed from clusters i and j, it is named min(i,j))
  mergelist=matrix(0,nblo-1,2)
  results=matrix(0,nblo-1,6)
  colnames(results)=c("merg1", "merg2", "new.clust", "agg.crit.hac","clust.crit.hac","S")
  idgr=matrix(0,nblo,nblo)
  idgr[nblo,]=1:nblo;
  ordr = c()
  fusions = -(1:nblo)
  hmerge =matrix(0,nblo-1,2)

  ncluster=nblo
  #for homogeneity
  quality=NULL

  for (level in 1:(nblo-1)) {

    ################
    # loops i et j
    ################

    # We search the two clusters that will be merged.

    deltamin=10^100
    for (i in 1:(ncluster-1)) {
      for (j in (i+1):(ncluster)) {
        # merging the i'th and j'th cluster:
        newcluster=c(which(cc==i), which(cc==j))
        Xj=list()
        for (tab in 1:length(newcluster)) {
          Xj[[tab]]=Xi[,,newcluster[tab]]
        }
        newcatatis=.crit_cataXj_fast(Xj, newcluster, S)
        deltacurrent=newcatatis-sum(crit[c(i,j)])
        # if deltacurrent is smaller than deltamin, the current information is saved:
        if (deltacurrent<deltamin) {
          deltamin=deltacurrent
          c1=cvar[cc==i]
          c2=cvar[cc==j]
          merge1=c(c1[1], c2[1])
          c1=code[cc==i]
          c2=code[cc==j]
          merge=c(c1[1], c2[1])
          cl_1=i
          cl_2=j
          catatismerge=newcatatis
        }
      }      # end of loop j
    }      # end of loop i
    ncluster=ncluster-1

    ###############################
    # renewal of the parameters
    ###############################

    Q_current=Q_current+deltamin
    mergelist[level,]=merge1
    results[level,1:5]=c(merge, nblo+level, deltamin, Q_current)
    crit[cl_1]=catatismerge
    crit=crit[-cl_2]
    cc[which(cc==cl_2)]=cl_1
    cc[which(cc>cl_2)]=cc[which(cc>cl_2)]-1
    cvar[cc==cl_1]=min(merge1)
    code[cc==cl_1]=nblo+level
    idgr[ncluster,]=cc
    indic= c(fusions[merge1[1]],fusions[merge1[2]])
    hmerge[level,]<-indic
    fusions[which(idgr[ncluster,]==cl_1)]<- level
    ordr <- .order_var(ordr,which(idgr[ncluster,]==cl_1))
    #overall homogeneity
    quality=c(quality, (nblo-Q_current)/nblo)

    if(printlevel==TRUE)
    {
      print(nblo-1-level)
    }
  }  # end of level


  # Dendogram:
  resultscah=list(labels=NameBlocks, height=results[,4], merge=hmerge,order=ordr)
  mytot<-resultscah
  class(mytot)="hclust"
  mydendC=as.dendrogram(mytot)

  #show the dendrogram
  if (Graph_dend==TRUE)
  {
    dev.new()
    cex=0.6
    par(cex=cex)
    par(mar=c(7, 4, 4, 2) + 0.1)
    plot(mydendC, type ="rectangle",  main="CLUSCATA Dendrogram", axes=TRUE, cex=cex,ylab="Height")
    par(cex=1)
    par(mar=c(5, 4, 4, 2) + 0.1)
  }


  #show the criterion and homogeneity evolutions
  if (Graph_bar==TRUE)
  {
    dev.new()
    barplot(results[,4][(nblo-gpmax):(nblo-1)],xlab="Nb clusters",ylab="delta", main="Variation of criterion",
            axisnames = TRUE,names.arg = paste((gpmax+1):2,"->",(gpmax):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

    dev.new()
    barplot(quality[(nblo-gpmax):(nblo-1)]*100,xlab="Nb clusters",ylab="Overall homogeneity (%)", main="Overall homogeneity (%)",
            axisnames = TRUE,names.arg = paste((gpmax):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")
  }


  #number of clusters advised
  criter=sort(results[,5],decreasing = TRUE)
  H=NULL
  for (k in 1:min(gpmax,nblo-2))
  {
    H[k]=(criter[k]/criter[k+1]-1)*(nblo-k-1)
  }
  nbgroup_hart=which.max(H[-(min(gpmax,nblo-2))]-H[-1])+1
  cat(paste("Recommended number of clusters =", nbgroup_hart),"\n")




  ######find the threshold rho for the noise cluster and consolidation
  res.consol=list()
  cutree_k=list()
  rho=NULL
  for (K in 1:gpmax)
  {
    coupe=cutree(mytot,K)
    cutree_k[[K]]=coupe
    if (Noise_cluster==TRUE)
    {
      if (is.null(rhoparam)==TRUE)
      {
        oldgroup=coupe
        #compute the compromises
        Ck=array(0,dim=c(nrow(Data),nvar,K))
        for (i in 1: K)
        {
          cluster=which(oldgroup==i)
          Xj=list()
          for (tab in 1:length(cluster)) {
            Xj[[tab]]=Xi[,,cluster[tab]]
          }
          catatisk=.crit_cataXj(Xj, cluster, S)
          Ck[,,i]=catatisk$C
        }

        #compute the criterions
        cr=list()
        cr2=list()
        for ( i in 1:nblo)
        {
          a=NULL
          for (k in 1:K)
          {
            C_k=as.matrix(Ck[,,k])
            normC=sum(diag(tcrossprod(C_k)))
            a=c(a,sum(diag(tcrossprod(Xi[,,i], C_k)))^2/(normC))
          }
          cr[[i]]=sqrt(a)
          if (K>1)
          {
            cr2[[i]]=sort(sqrt(a), decreasing = TRUE)[1:2]
          }else{
            cr2[[i]]=sqrt(a)
          }
        }
        rho[K]=mean(unlist(cr2))
      }else{
        rho[K]=rhoparam
      }
    } else{
      rho[K]=0
    }

    #consolidation

    res.consol[[K]]=cluscata_kmeans(Data, nblo, coupe, rho=rho[K], NameBlocks = NameBlocks, NameVar=NameVar, Itermax=Itermax,
                                    Graph_groups=FALSE, Warnings=Warnings)

  }
  names(cutree_k)=names(rho)=names(res.consol)=paste0("partition", 1:gpmax)

  #after consolidation
  diff_crit_bef=results[,4][(nblo-gpmax+1):(nblo-1)]

  if (Noise_cluster==FALSE)
  {
    #after consolidation
    overall_after=NULL
    crit_after=NULL
    for (i in 1: gpmax)
    {
      overall_after[i]=res.consol[[i]]$homogeneity[i+1,1]
      crit_after[i]=res.consol[[i]]$criterion
    }
    diff_crit_after=crit_after[-gpmax]-crit_after[-1]
    diff_crit_after=sort(diff_crit_after)
    overall_after=sort(overall_after, decreasing = TRUE)
    overall_homogeneity_ng=rbind(quality[(nblo-gpmax):(nblo-1)]*100, overall_after)
    diff_crit_ng=rbind(diff_crit_bef, diff_crit_after)
    rownames(diff_crit_ng)=rownames(overall_homogeneity_ng)=c("before consolidation", "after consolidation")
    colnames(diff_crit_ng)=paste((gpmax):2,"->",(gpmax-1):1)
    colnames(overall_homogeneity_ng)=paste (gpmax:1, "cluster(s)")
  }else{
    overall_homogeneity_ng=quality[(nblo-gpmax):(nblo-1)]*100
    diff_crit_ng=diff_crit_bef
    names(diff_crit_ng)=paste((gpmax):2,"->",(gpmax-1):1)
    names(overall_homogeneity_ng)=paste (gpmax:1, "cluster(s)")
  }



  #results
  res=c(res.consol, list(dend=mydendC, cutree_k=cutree_k,
                         overall_homogeneity_ng=round(overall_homogeneity_ng,1),
                         diff_crit_ng=round(diff_crit_ng,2),test_one_cluster=testonecluster, param=list(nblo=nblo, ng=nbgroup_hart,
                                                                                                        Noise_cluster=Noise_cluster, gpmax=gpmax, n=n, nvar=nvar), type="H+C"))

  class(res)="cluscata"

  return(res)

}
