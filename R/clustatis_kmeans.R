##=============================================================================

##' @title Compute the CLUSTATIS partitionning algorithm on different blocks of quantitative variables. Can be performed using a multi start strategy or initial partition provided by the user
##'
##'
##' @description
##' Partitionning of  of quantitative variables. Each cluster is associated with a compromise
##' computed by the STATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##' clustatis_kmeans(Data,Blocks, clust, nstart=40, rho=0, NameBlocks=NULL,
##' Itermax=30,Graph_groups=TRUE, Graph_weights=FALSE,
##'  scale=FALSE, print_attempt=FALSE)
##'
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data
##'
##' @param clust numerical vector or integer. Initial partition or number of starting partitions if integer. If numerical vector, the numbers must be 1,2,3,...,number of clusters
##'
##' @param nstart integer. Number of starting partitions. Default: 40
##'
##' @param rho numerical between 0 and 1. Threshold for the noise cluster. Default:0
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param Itermax numerical. Maximum of iterations by partitionning algorithm. Default: 30
##'
##' @param Graph_groups logical. Should each cluster compromise be plotted? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights in each cluster be plotted? Default: FALSE
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##' @param print_attempt logical. Print the number of remaining attempts in the multi-start case? Default: FALSE
##'
##'
##'
##'
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition. If Noise_cluster=TRUE, some blocks could be in the noise cluster ("K+1")
##'          \item rho: the threshold for the noise cluster
##'          \item homogeneity: percentage of homogeneity of the blocks in each cluster and the overall homogeneity
##'          \item rv_with_compromise: RV coefficient of each block with its cluster compromise
##'          \item weights: weight associated with each block in its cluster
##'          \item comp_RV: RV coefficient between the compromises associated with the various clusters
##'          \item compromise: the W compromise of each cluster
##'          \item coord: the coordinates of objects of each cluster
##'          \item inertia: percentage of total variance explained by each axis for each cluster
##'          \item rv_all_cluster: the RV coefficient between each block and each cluster compromise
##'          \item criterion: the CLUSTATIS criterion error
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##'
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods. Application to sensometrics. Food Quality and Preference, in Press.\cr
##' Llobell, F., Vigneau, E., Qannari, E. M. (2019). Clustering datasets by means of CLUSTATIS with identification of atypical datasets. Application to sensometrics. Food Quality and Preference, 75, 97-104.
##'
##'
##'
##' @keywords quantitative
##'
##' @examples
##'
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  #with multi-start
##'  cl_km=clustatis_kmeans(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks, clust=3)
##'  #with an initial partition
##'  cl=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks,
##'  Graph_dend=FALSE)
##'  partition=cl$cutree_k$partition3
##'  cl_km2=clustatis_kmeans(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks,
##'  clust=partition, Graph_weights=FALSE, Graph_groups=FALSE)
##'
##' @seealso   \code{\link{plot.clustatis}}, \code{\link{clustatis}}, \code{\link{summary.clustatis}}, \code{\link{statis}}
##'
##' @export


##=============================================================================








clustatis_kmeans=function(Data,Blocks, clust, nstart=40, rho=0, NameBlocks=NULL,Itermax=30,Graph_groups=TRUE,
                          Graph_weights=FALSE, scale=FALSE, print_attempt=FALSE){

  nblo=length(Blocks)

  #initialisation or muti-start?
  if(length(clust)==1)
  {
    ngroups=clust
    init=FALSE
  }else{
    ngroups=length(unique(clust))
    partition=clust
    init=TRUE
    if(sum(sort(unique(clust))!=(1:ngroups))>0)
    {
      stop("The clusters must be 1,2,...,number of clusters")
    }
  }

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

  #parapet for number of groups and number of blocks
  if(ngroups>nblo)
  {
    stop("number of groups > number of blocks")
  }

  #centering and scaling if necessary
  Data=scale(Data,center=TRUE,scale=scale)

  J=rep(1:nblo , times =  Blocks )      # indicates which block each variable belongs to

  #parapet for constant configurations and compute of Wi

  Wi=array(0,dim=c(nrow(Data),nrow(Data),nblo))
  for (i in 1:nblo)
  {
    Xi=as.matrix(Data[,J==i])
    Wi[,,i]=Xi%*%t(Xi)
    nor=sqrt(sum(diag(Wi[,,i]%*%Wi[,,i])))
    if(nor==0)
    {
      stop(paste("Configuration",i,"is constant"))
    }
    Wi[,,i]=Wi[,,i]/sqrt(sum(diag(Wi[,,i]%*%Wi[,,i])))#standardization
  }



  if(init==FALSE)
  {
    sumRV2=0
    for (tent in 1: nstart)
    {
      coupe=sample(1:ngroups,nblo,replace=TRUE)
      while(nlevels(as.factor(coupe))<ngroups)
      {
        coupe=sample(1:ngroups,nblo,replace=TRUE)#avoid empty groups
      }


      iter=0
      sumRV2_act=0
      goout=0
      oldgroup=coupe
      erreur=NULL
      while (goout!=1 & iter<=Itermax)
      {
        #avoid empty group
        if (nlevels(as.factor(oldgroup))<ngroups)
        {
          coupe=sample(1:ngroups,nblo,replace=TRUE)
          while(nlevels(as.factor(coupe))<ngroups)
          {
            coupe=sample(1:ngroups,nblo,replace=TRUE)
          }
          oldgroup=coupe
        }
        #compute the compromises
        Wk=array(0,dim=c(nrow(Data),nrow(Data),ngroups))
        lamb=NULL
        alpha=list()
        for (i in 1: ngroups)
        {
          cluster=which(oldgroup==i)
          Wj=list()
          for (tab in 1:length(cluster)) {
            Wj[[tab]]=Wi[,,cluster[tab]]
          }
          statisgr=.crit_statisWj(Wj)
          Wk[,,i]=statisgr$W
          lamb[i]=statisgr$lambda/length(cluster)
          alpha[[i]]=statisgr$u
          erreur[i]=statisgr$Q
        }


        #compute the criterion
        cr=list()
        for ( i in 1:nblo)
        {
          a=NULL
          for (k in 1: ngroups)
          {
            W_k=as.matrix(Wk[,,k])
            normW=sum(diag(W_k%*%W_k))
            W_i=Wi[,,i]
            a=c(a,sum(diag(W_i%*%W_k))^2/(normW))
          }
          cr[[i]]=sqrt(a)
        }

        #choice of the clusters
        newgroup=NULL
        for (i in 1:nblo)
        {
          if (max(cr[[i]])<rho)
          {
            newgroup[i]=0
          } else{
            newgroup[i]=which.max(cr[[i]])
          }
        }

        if (sum(newgroup!=oldgroup)>0)
        {
          oldgroup=newgroup
        }else{
          goout=1
        }

        #1 cluster left?
        if(0%in% oldgroup & nlevels(factor(oldgroup))!=ngroups+1)
        {
          warning(paste("One cluster is empty with "), ngroups," clusters, rho is too high")
          oldgroup[oldgroup==0]="K+1"
          return(oldgroup)
        }

        iter=iter+1
      }

      #best partition?
      quali_tent=sum(erreur)
      for (i in 1:nblo)
      {
        if (max(cr[[i]])<rho)
        {
          sumRV2_act=sumRV2_act+rho**2
        }else{
          sumRV2_act=sumRV2_act+max(cr[[i]]**2)
        }
      }
      if((sumRV2_act)>sumRV2)
      {
        sumRV2=sumRV2_act
        qual=quali_tent
        parti=oldgroup
        alpha_bon=alpha
        lamb_bon=lamb
        Wk_bon=Wk
        cr_bon=cr
        err_bon=sum(erreur)
        erreur_bon=erreur
      }
      if(print_attempt==TRUE)
      {
        print(nstart-tent)
      }
    }
  }

  if(init==TRUE)
  {
    sumRV2=0
    coupe=as.numeric(partition)
    iter=0
    sumRV2_act=0
    goout=0
    oldgroup=coupe
    erreur=NULL
    while (goout!=1 & iter<=Itermax)
    {
      #compute the compromises
      Wk=array(0,dim=c(nrow(Data),nrow(Data),ngroups))
      lamb=NULL
      alpha=list()
      for (i in 1: ngroups)
      {
        cluster=which(oldgroup==i)
        Wj=list()
        for (tab in 1:length(cluster)) {
          Wj[[tab]]=Wi[,,cluster[tab]]
        }
        statisgr=.crit_statisWj(Wj)
        Wk[,,i]=statisgr$W
        lamb[i]=statisgr$lambda/length(cluster)
        alpha[[i]]=statisgr$u
        erreur[i]=statisgr$Q
      }


      #compute the criterion
      cr=list()

      for ( i in 1:nblo)
      {
        a=NULL
        for (k in 1: ngroups)
        {
          W_k=as.matrix(Wk[,,k])
          normW=sum(diag(W_k%*%W_k))
          W_i=Wi[,,i]
          a=c(a,sum(diag(W_i%*%W_k))^2/(normW))
        }
        cr[[i]]=sqrt(a)
      }

      #choice of the clusters
      newgroup=NULL
      for (i in 1:nblo)
      {
        if (max(cr[[i]])<rho)
        {
          newgroup[i]=0
        }else{
          newgroup[i]=which.max(cr[[i]])
        }
      }

      if (sum(newgroup!=oldgroup)>0)
      {
        oldgroup=newgroup
      }else{
        goout=1
      }

      #1 cluster left?
      if(0%in% oldgroup & nlevels(factor(oldgroup))!=ngroups+1)
      {
        warning(paste("One cluster is empty with "), ngroups," clusters, rho is too high")
        oldgroup[oldgroup==0]="K+1"
        return(oldgroup)
      }
      iter=iter+1
    }

    #results of the partition
    quali_tent=sum(erreur)
    for (i in 1:nblo)
    {
      if (max(cr[[i]])<rho)
      {
        sumRV2_act=sumRV2_act+rho**2
      }else{
        sumRV2_act=sumRV2_act+max(cr[[i]]**2)
      }
    }

    sumRV2=sumRV2_act
    qual=quali_tent
    parti=oldgroup
    alpha_bon=alpha
    lamb_bon=lamb
    Wk_bon=Wk
    cr_bon=cr
    err_bon=sum(erreur)
    erreur_bon=erreur
  }

  #to simplify names
  oldgroup=parti
  Wk=Wk_bon
  alpha=alpha_bon
  lamb=lamb_bon
  cr=cr_bon

  #tab1: homogeneity
  Wj=list()
  for (i in 1:nblo)
  {
    Wj[[i]]=Wi[,,i]
  }
  stati=.crit_statisWj(Wj)
  lambda_tot=stati$lambda/nblo
  ne=NULL
  for (i in 1:ngroups)
  {
    ne=c(ne,length(which(oldgroup==i)))
  }
  adios=which(oldgroup==0)
  overall=(sum(lamb*ne))/(sum(ne))
  if(length(adios)==0)
  {
    tab1=data.frame(homogeneity=round(c(lamb, overall, lambda_tot),3)*100,nb_blocks=c(ne, nblo, nblo))
    rownames(tab1)=c(paste("Cluster",1:ngroups),"Overall","One group")
  }else{
    cluster=which(oldgroup==0)
    Wj=list()
    for (tab in 1:length(cluster)) {
      Wj[[tab]]=Wi[,,cluster[tab]]
    }
    statisk1=.crit_statisWj(Wj)
    lambk1=statisk1$lambda/length(cluster)
    tab1=data.frame(homogeneity=round(c(lamb,lambk1, overall,lambda_tot),3)*100, nb_blocks=c(ne,length(adios), sum(ne), nblo))
    rownames(tab1)=c(paste("Cluster",1:ngroups), "Noise cluster", "Overall", "One group")
  }
colnames(tab1)[1]="homogeneity (%)"

#rv with compromise and weights
  rvcomp=list()
  for (k in 1:ngroups)
  {
    cluster=which(parti==k)
    temp=NULL
    for (i in cluster)
    {
      temp=c(temp, round(max(cr[[i]]),3))
    }
    rvcomp[[k]]=temp
    alpha[[k]]=as.vector(alpha[[k]])
    names(rvcomp[[k]])=names(alpha[[k]])=NameBlocks[cluster]
  }

  #groups
  partit=parti
  parti[adios]="K+1"
  names(parti)=NameBlocks
  grou=data.frame(Cluster=parti)

  #RV among compromises
  tab4=matrix(0,ngroups,ngroups)
  for (i in 1: ngroups)
  {
    for(j in 1:ngroups)
    {
      W_k=as.matrix(Wk_bon[,,i])
      W_l=as.matrix(Wk_bon[,,j])
      normWk=sqrt(sum(diag(W_k%*%W_k)))
      normWl=sqrt(sum(diag(W_l%*%W_l)))
      tab4[i,j]=sum(diag(W_l%*%W_k))/(normWl*normWk)
    }
  }

  rownames(tab4)=paste("Cluster",1:ngroups)
  colnames(tab4)=paste("Cluster",1:ngroups)

  #proximity inter-cluster
  prox=eigen(tab4)$values[1]/ngroups

  #coordinates
  compromise=list()
  coord=list()
  inertia=list()
  for (i in 1:ngroups)
  {
    e=svd(Wk[,,i])
    C=e$u%*%sqrt(diag(abs(e$d)))
    C=C[,-ncol(C)]
    if(i==1)
    {
      Cref=C
    }
    if (i!=1)
    {
      for (j in 1:ncol(C))
      {
        if(var(C[,j])>10^(-12) & var(Cref[,j])>10^(-12))
        {
          if (cor(Cref[,j],C[,j])<0)
          {
            C[,j]=-C[,j]#axes orientation
          }
        }
      }
    }
    pouriner=round(e$d/sum(e$d)*100,2)
    comprom=Wk[,,i]
    rownames(comprom)=colnames(comprom)=rownames(C)=rownames(Data)
    compromise[[i]]=comprom
    colnames(C)=paste("Dim", 1:(nrow(Data)-1))
    coord[[i]]=C
    pouriner=pouriner[-length(pouriner)]
    names(pouriner)=paste("Dim", 1:(nrow(Data)-1))
    inertia[[i]]=pouriner

  #graphical representations
  if (Graph_groups==TRUE)
  {
      dev.new()
      plot(C[,1],C[,2],type="n",lwd=5,pch=16,asp=1,xlab=paste("Dim 1 (",pouriner[1],"%)"), ylab=paste("Dim 2 (",pouriner[2],"%)"),xlim=c(min(C[,1])-0.05,max(C[,1])+0.05))
      text(C[,1],C[,2],rownames(Data),col=rainbow(nrow(Data)))
      abline(h=0,v=0)
      title(paste("Cluster",i))
    }
  }
  names(compromise)=names(rvcomp)=names(alpha)=paste("Cluster",1:ngroups)

  if(length(adios)==0)
  {
    criterion=err_bon
  }else{
    criterion=NULL
  }

  #weights barplot
  if(Graph_weights==TRUE)
  {
    for (i in 1:ngroups)
    {
      dev.new()
      barplot(alpha[[i]])
      title(paste("Weights in the cluster",i))
    }
  }

  #rounds
  for (i in 1: ngroups)
  {
    compromise[[i]]=round(compromise[[i]],3)
    coord[[i]]=round(coord[[i]],3)
  }
  for (i in 1:nblo)
  {
    cr[[i]]=round(cr[[i]],3)
    Wi[,,i]=round(Wi[,,i],3)
    names(cr[[i]])=paste("Cluster", 1:ngroups)
  }
  if(nblo>1)
  {
    dimnames(Wi)[[1]]=dimnames(Wi)[[2]]=rownames(Data)
    dimnames(Wi)[[3]]=NameBlocks
  }




  names(cr)=NameBlocks
  names(compromise)=names(coord)=names(inertia)=paste("Cluster",1:ngroups)



  #results
  res=list(group=grou, rho=rho, homogeneity=tab1, rv_with_compromise=rvcomp, weights=alpha, compRV=round(tab4,2),
           compromise=compromise, coord=coord, inertia=inertia,
           rv_all_cluster=cr, criterion=criterion, param=list(nblo=nblo, ng=ngroups, n=n), type="K")
  class(res)="clustatis"

  return(res)


}
