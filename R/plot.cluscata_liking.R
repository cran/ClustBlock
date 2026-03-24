## =============================================================================

##' @title Displays the CLUSCATA-liking graphs
##'
##' @usage
##' \method{plot}{cluscata_liking}(x, DataCata, Dataliking, ngroups=NULL,
##'  center=TRUE, scale=FALSE,
##'  sep_graphs=TRUE, Graph_dend=TRUE, Graph_bar=FALSE,
##'   cex=1,  xlimfreq=c(0,0.6), ylimchangelik=c(-2, 2), ...)
##' @description
##' This function plots dendrogram, variation of the merging criterion, weights and CATATIS map of each cluster
##'
##' @param x object of class 'cluscata_liking'.
##'
##' @param DataCata data frame or matrix where the blocks of binary variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param Dataliking data frame or matrix where the products are in rows and the assessors in columns
##'
##' @param ngroups number of groups to consider. Default: recommended  number of clusters
##'
##' @param center Centering of consumer liking. Default: TRUE
##'
##' @param scale Scaling  of consumer liking. Default: FALSE
##'
##' @param sep_graphs logical. Should cata and liking data analyzed separately first?
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Also available after consolidation if Noise_cluster=FALSE. Default: FALSE
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param xlimfreq vector of 2 values. Used for the graph of the frequency of citations
##'
##' @param ylimchangelik vector of 2 values. Used for the graphs of change in liking
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return the CLUSCATA-liking graphs
##'
##'
##' @keywords CATA liking
##'
##' @examples
##' \donttest{
##' data(cata_ryebread)
##' data(liking_ryebread)
##' cataliking=combinCATALiking(cata_ryebread, liking_ryebread)
##'
##'#with only 40 subjects
##' resclustcatal=cluscata_liking(Data=cataliking[,1:(40*14)], nblo=40, gpmax=5)
##' plot(resclustcatal, cata_ryebread[,1:(40*14)], liking_ryebread[,1:40],
##' xlimfreq=c(0,0.6),  ylimchangelik=c(-2, 2))
##'
##' }
##'
##' @seealso   \code{\link{cluscata_liking}} , \code{\link{summary.cluscata_liking}}
##'
##' @export


## =============================================================================


plot.cluscata_liking=function(x, DataCata, Dataliking, ngroups=NULL, center=TRUE,
                              scale=FALSE, sep_graphs=TRUE, Graph_dend=TRUE,
                              Graph_bar=FALSE, cex=1, xlimfreq=c(0,0.6),
                              ylimchangelik=c(-2, 2),...)
{
  res.cluscata=x

  if (inherits(res.cluscata, "cluscata_liking") == FALSE) {
    stop("The class of the object must be 'cluscata_liking'")
  }

  liking=as.matrix(Dataliking)
  DataCata=as.matrix(DataCata)
  n=res.cluscata$param$n

  if(is.null(ngroups) | res.cluscata$type=="K")
  {
    ngroups=res.cluscata$param$ng
  }

  if(res.cluscata$type=="H+C")
  {
    if(ngroups>res.cluscata$param$gpmax)
    {
      stop("ngroups>gpmax")
    }
  }

  v_colors <- c("blue","red","green","black",
                "purple","orange","yellow","tomato","pink",
                "gold","cyan","turquoise","violet",
                "khaki","maroon","chocolate","burlywood")

  #show the dendrogram
  if (Graph_dend==TRUE & res.cluscata$type=="H+C")
  {
    dev.new()
    par(cex=0.6)
    par(mar=c(7, 4, 4, 2) + 0.1)
    plot(res.cluscata$dend, type ="rectangle",  main="CLUSCATA-liking Dendrogram", axes=TRUE, cex=cex,ylab="Height")
    par(cex=1)
    par(mar=c(5, 4, 4, 2) + 0.1)
  }



  #show the criterion and homogeneity evolutions
  if (Graph_bar==TRUE & res.cluscata$type=="H+C")
  {
    nblo=res.cluscata$param$nblo
    gpmax=res.cluscata$param$gpmax

    dif=res.cluscata$diff_crit_ng[1,]
    quality=res.cluscata$overall_homogeneity_ng[1,]


    dev.new()
    barplot(dif,xlab="Nb clusters",ylab="delta", main="Variation of criterion before consolidation",
            axisnames = TRUE,names.arg = paste((gpmax):2,"->",(gpmax-1):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

       dif=res.cluscata$diff_crit_ng[2,]
    quality=res.cluscata$overall_homogeneity_ng[2,]

    dev.new()
    barplot(dif,xlab="Nb clusters",ylab="delta", main="Variation of criterion after consolidation",
            axisnames = TRUE,names.arg = paste((gpmax):2,"->",(gpmax-1):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

    }

  if(res.cluscata$type=="H+C")
  {
    res.cluscata=res.cluscata[[ngroups]]
  }

  nblo=res.cluscata$param$nblo
  nvar=res.cluscata$param$nvar
  nprod=res.cluscata$param$n
  Blocks=rep(nvar,nblo)
  J=rep(1:nblo , times =  Blocks )
  likingc=scale(liking, center = center, scale=scale)

  #show analyses of CATA data and liking data separately

  if(sep_graphs==TRUE)
  {

    #limits
    ylminref=0
    ylmaxref=0
    for (k in 1:ngroups)
    {
      #liking
      cluster=which(res.cluscata$group==k)
      likingCluster=likingc[,cluster]
      likingMoy=apply(likingCluster, 1, "mean")
      ylmin=min(likingMoy)
      if (ylmin < ylminref)
      {
        ylminref = ylmin
      }
      ylmax=max(likingMoy)
      if (ylmax > ylmaxref)
      {
        ylmaxref = ylmax
      }
    }

    for (k in 1:ngroups)
    {
      #cata
      cluster=which(res.cluscata$group==k)
      DataCataCluster=DataCata[,J%in%cluster]
      res.cat=catatis(DataCataCluster, length(cluster), Graph_weights = FALSE, Graph = FALSE)
      plot(res.cat, tit=paste("CATATIS analysis of Cluster", k), Graph_weights = FALSE, col.obj=v_colors[1:res.cluscata$param$n], col.attr="blue4")

      #liking
      likingCluster=likingc[,cluster]
      likingMoy=apply(likingCluster, 1, "mean")
      ylmin=min(likingMoy)
      ylmax=max(likingMoy)
      dev.new()
      barplot(likingMoy, col=v_colors[1:res.cluscata$param$nvar],
              main=paste("Means of liking in Cluster", k),
              ylim=c(ylminref, ylmaxref))

    }
  }

  #show analyses of CATA data and liking data together

  for (k in 1:ngroups)
  {
    cluster=which(res.cluscata$group==k)
    DataCataCluster=DataCata[,J%in%cluster]
    likingCluster=likingc[,cluster]
    Diff=NULL
    for (i in 1:nvar)
    {
      DataAttr=NULL
      for (j in seq(i,nvar*length(cluster), nvar))
      {
        DataAttr=c(DataAttr, DataCataCluster[,j])
      }
      Resume=cbind(DataAttr, as.vector(likingCluster))

      moyCheck=mean(Resume[which(Resume[,1]==1),2])
      moyNotCheck=mean(Resume[which(Resume[,1]==0),2])
      Diff=c(Diff,moyCheck-moyNotCheck)

    }
    names(Diff)=colnames(DataCata)[1:nvar]

    dev.new()
    barplot(Diff, col=v_colors[1:res.cluscata$param$nvar],
            main=paste("Change in liking in cluster", k), ylim=ylimchangelik)

    dev.new()
    cube=array(0, dim=c(nprod, nvar, length(cluster)))
    cpt=0
    for (i in cluster)
    {
      cpt=cpt+1
      cube[,,cpt]=DataCata[,J==i]
    }
    freqcit=apply(cube, 2, mean)
    plot(freqcit, Diff, type="n", xlab="freq of citation",
         ylab="Change in liking", main = paste("Cluster", k),
         xlim=xlimfreq, ylim=ylimchangelik)
    text(freqcit, Diff, dimnames(DataCata)[[2]][1:nvar])
    abline(h=0)

  }



}
