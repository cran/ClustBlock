##=============================================================================


##' @title Compute the CLUSCATA partitioning algorithm on different blocks from a CATA experiment. Can be performed using a multi-start strategy or initial partition provided by the user.
##'
##'
##' @description
##' Partitioning of binary Blocks from a CATA experiment. Each cluster is associated with a compromise computed by the CATATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##'cluscata_kmeans(Data,nblo, clust, nstart=100, rho=0, NameBlocks=NULL, NameVar=NULL,
##'                Itermax=30, Graph_groups=TRUE, print_attempt=FALSE, Warnings=FALSE)
##'
##'
##' @param Data data frame or matrix where the blocks of binary variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param nblo  numerical. Number of blocks (subjects).
##'
##' @param clust numerical vector or integer. Initial partition or number of starting partitions if integer. If numerical vector, the numbers must be 1,2,3,...,number of clusters
##'
##' @param nstart numerical. Number of starting partitions. Default: 100
##'
##' @param rho numerical between 0 and 1. Threshold for the noise cluster. If 0, there is no noise cluster. Default: 0
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the number of blocks. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param NameVar string vector. Name of each variable (attribute, the same names for each subject). Length must be equal to the number of attributes. If NULL, the colnames of the first block are taken. Default: NULL
##'
##' @param Itermax numerical. Maximum of iterations by partitionning algorithm. Default: 30
##'
##' @param Graph_groups logical. Should each cluster compromise graphical representation be plotted? Default: TRUE
##'
##' @param print_attempt logical. Print the number of remaining attempts in multi-start case? Default: FALSE
##'
##' @param Warnings logical. Display warnings about the fact that none of the subjects in some clusters checked an attribute or product? Default: FALSE
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition. If rho>0, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold for the noise cluster
##'          \item homogeneity: percentage of homogeneity of the subjects in each cluster and the overall homogeneity
##'          \item s_with_compromise: Similarity coefficient of each subject with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item compromise: The compromise of each cluster
##'          \item CA: The correspondance analysis results on each cluster compromise (coordinates, contributions...)
##'          \item inertia: percentage of total variance explained by each axis of the CA for each cluster
##'          \item s_all_cluster: the similarity coefficient between each subject and each cluster compromise
##'          \item param: parameters called
##'          \item criterion: the CLUSCATA criterion error
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
##' cl_km=cluscata_kmeans(Data=straw[,1:(16*40)], nblo=40, clust=3)
##' #plot(cl_km, Graph_groups=FALSE, Graph_weights = TRUE)
##' summary(cl_km)
##'}
##'
##' @seealso   \code{\link{plot.cluscata}} , \code{\link{summary.cluscata}}, \code{\link{catatis}}, \code{\link{cluscata}}, \code{\link{change_cata_format}}
##'
##' @export

##=============================================================================



cluscata_kmeans=function(Data,nblo, clust, nstart=100, rho=0, NameBlocks=NULL, NameVar=NULL, Itermax=30,
                         Graph_groups=TRUE, print_attempt=FALSE, Warnings=FALSE){



  #initialisation or muti-start?
  if(length(clust)==1)
  {
    ngroups=clust
    init=FALSE
  }else{
    ngroups=length(unique(clust))
    partition=clust
    init=TRUE
  }

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
  if(is.null(colnames(Data))) colnames(Data)=paste0("Y",1:Blocks[1])

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
    stop("At least 3 objects are required")
  }

  #parapet for number of blocks
  if(nblo<4)
  {
    stop("At least 4 blocks are required")
  }

  #parapet for number of groups and number of blocks
  if(ngroups>nblo)
  {
    stop("number of groups > number of blocks")
  }


  #no NA
  if(sum(is.na(Data))>0)
  {
    print("NA detected:")
    tabna=which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }


  #compute of Xi
  Xi=array(0,dim=c(nrow(X),nvar,nblo))
  for (i in 1:nblo)
  {
    Ai=as.matrix(X[,J==i])
    nor=sqrt(sum(diag(tcrossprod(Ai,Ai))))
    if(nor==0)
    {
      stop(paste("error: the subject",NameBlocks[i], "has only 0"))
    }
    Xi[,,i]=Ai/nor
  }


  # S matrix:
  S=matrix(0,nblo,nblo)
  diag(S)=rep(1,nblo)
  for (i in 1:(nblo-1)) {
    for (j in (i+1):nblo) {
      S[i,j]=sum(diag(tcrossprod(Xi[,,i],Xi[,,j])))
      S[j,i]=S[i,j]
    } }

  if(init==FALSE)
  {
    qual=10000
    sums2=0
    for (tent in 1: nstart)
    {
      coupe=sample(1:ngroups,nblo,replace=TRUE)
      while(nlevels(as.factor(coupe))<ngroups)
      {
        coupe=sample(1:ngroups,nblo,replace=TRUE)#avoid empty groups
      }


      iter=0
      sums2_act=0

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
        #compute of the compromises
        Ck=array(0,dim=c(nrow(X),nvar,ngroups))
        lamb=NULL
        alpha=list()
        for (i in 1: ngroups)
        {
          cluster=which(oldgroup==i)
          Xj=list()
          for (tab in 1:length(cluster)) {
            Xj[[tab]]=Xi[,,cluster[tab]]
          }
          catatisgr=.crit_cataXj(Xj, cluster, S)
          Ck[,,i]=catatisgr$C
          lamb[i]=catatisgr$lambda/length(cluster)
          alpha[[i]]=catatisgr$alpha
          erreur[i]=catatisgr$Q
        }

        #compute the criterion
        cr=list()
        for ( i in 1:nblo)
        {
          a=NULL
          for (k in 1: ngroups)
          {
            C_k=as.matrix(Ck[,,k])
            normC=sum(diag(tcrossprod(C_k)))
            X_i=Xi[,,i]
            a=c(a,sum(diag(tcrossprod(X_i, C_k)))^2/(normC))
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
        iter=iter+1
        #1 cluster left?
        if(0%in% oldgroup & nlevels(factor(oldgroup))!=ngroups+1)
        {
          warning(paste("One cluster is empty with "), ngroups," clusters, rho is too high")
          oldgroup[oldgroup==0]="K+1"
          return(oldgroup)
        }
      }

      #best partition?
      quali_tent=sum(erreur)
      for (i in 1:nblo)
      {
        if (max(cr[[i]])<rho)
        {
          sums2_act=sums2_act+rho**2
        }
        else
        {
          sums2_act=sums2_act+max(cr[[i]]**2)
        }
      }
      if((sums2_act)>sums2)
      {
        sums2=sums2_act
        qual=quali_tent
        parti=oldgroup
        alpha_bon=alpha
        lamb_bon=lamb
        Ck_bon=Ck
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
    sums2=0
    qual=100000
    coupe=as.numeric(partition)
    iter=0
    sums2_act=0
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
      #compute of the compromises
      Ck=array(0,dim=c(nrow(X),nvar,ngroups))
      lamb=NULL
      erreur=NULL
      alpha=list()
      for (i in 1: ngroups)
      {
        cluster=which(oldgroup==i)
        Xj=list()
        for (tab in 1:length(cluster)) {
          Xj[[tab]]=Xi[,,cluster[tab]]
        }
        catatisgr=.crit_cataXj(Xj, cluster, S)
        Ck[,,i]=catatisgr$C
        lamb[i]=catatisgr$lambda/sum(diag(catatisgr$S))
        alpha[[i]]=catatisgr$alpha
        erreur[i]=catatisgr$Q
      }


      #compute the criterion
      cr=list()
      for ( i in 1:nblo)
      {
        a=NULL
        for (k in 1: ngroups)
        {
          C_k=as.matrix(Ck[,,k])
          normC=sum(diag(tcrossprod(C_k)))
          X_i=Xi[,,i]
          a=c(a,sum(diag(tcrossprod(X_i, C_k)))^2/(normC))
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
        }
        else
        {
          newgroup[i]=which.max(cr[[i]])
        }
      }



      if (sum(newgroup!=oldgroup)>0)
      {
        oldgroup=newgroup
      }
      else
      {
        goout=1
      }
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
        sums2_act=sums2_act+rho**2
      }else{
        sums2_act=sums2_act+max(cr[[i]]**2)
      }
      sums2=sums2_act
      qual=quali_tent
      parti=oldgroup
      alpha_bon=alpha
      lamb_bon=lamb
      Ck_bon=Ck
      cr_bon=cr
      err_bon=sum(erreur)
      erreur_bon=erreur
    }
  }


  #to simplify names
  oldgroup=parti
  Ck=Ck_bon
  alpha=alpha_bon
  lamb=lamb_bon
  cr=cr_bon


  #tab1: homogeneity
  Xj=list()
  for (i in 1:nblo)
  {
    Xj[[i]]=Xi[,,i]
  }
  catati=.crit_cataXj(Xj, 1:nblo, S)
  lambda_tot=catati$lambda/nblo
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
    Xj=list()
    for (tab in 1:length(cluster)) {
      Xj[[tab]]=Xi[,,cluster[tab]]
    }
    catatisk1=.crit_cataXj(Xj, cluster, S)
    lambk1=catatisk1$lambda/length(cluster)
    tab1=data.frame(homogeneity=round(c(lamb,lambk1, overall,lambda_tot),3)*100, nb_blocks=c(ne,length(adios), sum(ne), nblo))
    rownames(tab1)=c(paste("Cluster",1:ngroups), "Noise cluster", "Overall", "One group")
  }
  colnames(tab1)[1]="homogeneity (%)"


  #s with compromise and weights
  scomp=list()
  for (k in 1:ngroups)
  {
    cluster=which(parti==k)
    temp=NULL
    for (i in cluster)
    {
      temp=c(temp, round(max(cr[[i]]),3))
    }
    scomp[[k]]=temp
    alpha[[k]]=as.vector(alpha[[k]])
    names(scomp[[k]])=names(alpha[[k]])=NameBlocks[cluster]
  }
  #groups
  partit=parti
  parti[adios]="K+1"
  names(parti)=NameBlocks
  grou=data.frame(Cluster=parti)

  #compute of the coordinates

  compro=list()
  AFC=list()
  inertia=list()
  eigenvalues=list()
  for (i in 1:ngroups)
  {
    compromis=Ck[,,i]
    rownames(compromis)=rownames(X)
    if (is.null(NameVar)==TRUE)
    {
      colnames(compromis)=colnames(X)[1:nvar]
    }else{
      colnames(compromis)=NameVar
    }
    compro[[i]]=compromis
    #CA
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
      if (Warnings==TRUE)
      {
        warning("Partion in ", paste(ngroups), " clusters: no subject in cluster ", paste(i), " checked the attribute(s): ", paste(colnames(compromis)[colomnnull], collapse=","))
      }
    }
    if(length(rownull)>0)
    {
      compromis2=compromis2[-rownull,]
      if (Warnings==TRUE)
      {
        warning("Partion in ", paste(ngroups), " clusters: no subject in cluster ", paste(i), " has a 1 for the product(s):  ", paste(rownames(compromis)[rownull], collapse=","))
      }
    }

    e=CA(compromis2,graph=FALSE)

    #inertia of axes
    pouriner=round(e$eig[,2],2)
    eigenvalues=round(e$eig[,1],4)
    inertia[[i]]=pouriner

    C=e$row$coord
    if(i==1)
    {
      if(length(rownull)==0)
      {
        Cref=C
        pb=0
      }else{
        pb=1
      }
    }
    if (i!=1 & length(rownull)==0 & pb==0)
    {
      for (j in 1:ncol(C))
      {
        if(var(C[,j])>10^(-12) & var(Cref[,j])>10^(-12))
        {
          if (cor(Cref[,j],C[,j])<0)
          {
            e$row$coord[,j]=-e$row$coord[,j]#axes orientation
            e$col$coord[,j]=-e$col$coord[,j]
          }
        }
      }
    }

    AFC[[i]]=e

    #graphical representation of each cluster
    if (Graph_groups==TRUE)
    {

      dev.new()
      options(ggrepel.max.overlaps = Inf)
      print(plot.CA(e,title=paste("Cluster",i)))
    }
  }



  compromise=compro
  for (i in 1: ngroups)
  {
    compromise[[i]]=round(compromise[[i]],2)
  }
  names(compromise)=paste("Cluster",1:ngroups)

  if(length(adios)==0)
  {
    criterion=err_bon
  }else{
    criterion=NULL
  }

  for (i in 1:nblo)
  {
    cr[[i]]=round(cr[[i]],3)
    names(cr[[i]])=paste("Cluster", 1:ngroups)
  }


  names(AFC)=names(scomp)=names(alpha)=names(inertia)=paste("Cluster",1:ngroups)
  names(cr)=NameBlocks

  #results
  res=list(group=grou,rho=rho, homogeneity=tab1, s_with_compromise=scomp, weights=alpha,
           compromise=compromise, CA=AFC, inertia=inertia,
           s_all_cluster=cr, criterion=criterion, param=list(nblo=nblo, ng=ngroups, n=n, nvar=nvar), type="K")

  class(res)="cluscata"

  return(res)

}
