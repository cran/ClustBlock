.cluscata_kmeans_liking=function(Data, nblo, clust,  NameBlocks=NULL, NameVar=NULL,
                                Itermax=30, print_attempt=FALSE){

  ngroups=length(unique(clust))
  partition=clust

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
  Xi=array(0,dim=c(nrow(X),Blocks[1],nblo))
  for (i in 1:nblo)
  {
    Ai=as.matrix(X[,J==i])
    nor=sqrt(sum(diag(Ai%*%t(Ai))))
    if(nor==0)
    {
      stop(paste("error: the subject",NameBlocks[i], "has only 0"))
    }
    Xi[,,i]=Ai/nor
  }

  ###partitioning algorithm
  sums2=0
  coupe=as.numeric(partition)
  iter=0
  sums2_act=0
  goout=0
  oldgroup=coupe
  erreur=NULL

  while (goout!=1 & iter<=Itermax)
  {
    #compute of the compromises
    Ck=array(0,dim=c(nrow(X),Blocks[1],ngroups))
    erreur=NULL
    for (i in 1: ngroups)
    {
      cluster=which(oldgroup==i)
      Xj=list()
      for (tab in 1:length(cluster)) {
        Xj[[tab]]=Xi[,,cluster[tab]]
      }
      catatisgr=.crit_cataMean(Xj)
      Ck[,,i]=catatisgr$C
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
        X_i=Xi[,,i]
        er=X_i-C_k
        a=c(a,sum(diag(er%*%t(er))))
      }
      cr[[i]]=a
    }

    #choice of the clusters
    newgroup=NULL
    for (i in 1:nblo)
    {
      newgroup[i]=which.min(cr[[i]])
    }


    if (sum(newgroup!=oldgroup)>0)
    {
      oldgroup=newgroup
    }
    else
    {
      goout=1
    }
    iter=iter+1
  }

  for (i in 1:nblo)
  {
    parti=oldgroup
    Ck_bon=Ck
    cr_bon=cr
    err_bon=sum(erreur)
    erreur_bon=erreur
  }


  #to simplify names
  oldgroup=parti
  Ck=Ck_bon
  cr=cr_bon

  #groups
  partit=parti
  names(parti)=NameBlocks

  grou=data.frame(Cluster=parti)

  #compromises
  compro=list()
  for (i in 1:ngroups)
  {
    compromis=Ck[,,i]
    rownames(compromis)=rownames(X)
    if (is.null(NameVar)==TRUE)
    {
      colnames(compromis)=colnames(X)[1:(Blocks[1])]
    }else{
      colnames(compromis)=NameVar
    }
    compro[[i]]=compromis
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
      warning("Partion in ", paste(ngroups), " clusters: no subject in cluster ", paste(i), " checked the attribute(s): ", paste(colnames(compromis)[colomnnull], collapse=","))
    }
    if(length(rownull)>0)
    {
      compromis2=compromis2[-rownull,]
      warning("Partion in ", paste(ngroups), " clusters: no subject in cluster ", paste(i), " has a 1 for the product(s):  ", paste(rownames(compromis)[rownull], collapse=","))
    }

  }

  compromise=compro
  names(compromise)=paste("Cluster",1:ngroups)

  criterion=err_bon

  for (i in 1:nblo)
  {
    cr[[i]]=round(cr[[i]],3)
    names(cr[[i]])=paste("Cluster", 1:ngroups)
  }


  names(cr)=NameBlocks


  #results
  res=list(group=grou,
           compromise=compromise,
           dist_all_cluster=cr, criterion=criterion,
           param=list(nblo=nblo, ng=ngroups, n=n, nvar=nvar), type="K")

  class(res)="cluscata_liking"

  return(res)

}
