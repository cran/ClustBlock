##=============================================================================

##' @title Displays the CLUSTATIS graphs
##'
##' @usage
##' \method{plot}{clustatis}(x, ngroups=NULL, Graph_groups=TRUE, Graph_dend=TRUE,
##' Graph_bar=FALSE, Graph_weights=FALSE, axes=c(1,2), col=NULL, cex=1, font=1, ...)
##'
##' @description
##' This function plots dendrogram, variation of the merging criterion, weights and STATIS map of each cluster
##'
##'
##' @param x object of class 'clustatis'.
##'
##' @param ngroups number of groups to consider.  Ignored for clustatis_kmeans results. Default: recommended  number of clusters
##'
##' @param Graph_groups logical. Should each cluster compromise graphical representation be plotted? Default: TRUE
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Also available after consolidation if Noise_cluster=FALSE. Default: FALSE
##'
##' @param Graph_weights logical. Should the barplot of the weights in each cluster be plotted? Default: FALSE
##'
##' @param axes  numerical vector (length 2). Axes to be plotted. Default: c(1,2)
##'
##' @param col vector. Color for each object. Default: rainbow(nrow(Data))
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param font numerical. Integer specifying font to use for text. 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol. Default: 1
##'
##' @param ... further arguments passed to or from other methods
##
##'
##'
##' @return the CLUSTATIS graphs
##'
##'
##'
##' @keywords quantitative
##'
##' @examples
##' \donttest{
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  cl=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks)
##'  plot(cl, ngroups=3, Graph_dend=FALSE)
##'  plot(cl, ngroups=3,  Graph_dend=FALSE, axes=c(1,3))
##'  graphics.off()
##'  }
##'
##' @seealso   \code{\link{clustatis}} , \code{\link{clustatis_kmeans}}
##'
##' @export


##=============================================================================





plot.clustatis=function(x, ngroups=NULL, Graph_groups=TRUE, Graph_dend=TRUE,
                        Graph_bar=FALSE, Graph_weights=FALSE, axes=c(1,2), col=NULL,
                        cex=1,font=1, ...)
{
  res.clustatis=x
  if(class(res.clustatis)!="clustatis")
  {
    stop("The class of the object must be 'clustatis'")
  }

  n=res.clustatis$param$n

  if(is.null(col)==TRUE)
  {
    col=rainbow(n)
  }
  if(length(col)!= n)
  {
    stop("col must be the size of the number of objects")
  }

  if(is.null(ngroups) | res.clustatis$type=="K")
  {
    ngroups=res.clustatis$param$ng
  }

  if(res.clustatis$type=="H+C")
  {
     if(ngroups>res.clustatis$param$gpmax)
      {
        stop("ngroups>gpmax")
      }
  }

  #show the dendrogram
  if (Graph_dend==TRUE & res.clustatis$type=="H+C")
  {
    dev.new()
    par(cex=0.6)
    par(mar=c(7, 4, 4, 2) + 0.1)
    plot(res.clustatis$dend, type ="rectangle",  main="CLUSTATIS Dendrogram", axes=TRUE, cex=cex,ylab="Height")
    par(cex=1)
    par(mar=c(5, 4, 4, 2) + 0.1)
  }



  #show the criterion and homogeneity evolutions
  if (Graph_bar==TRUE & res.clustatis$type=="H+C")
  {
    nblo=res.clustatis$param$nblo
    Noise_cluster=res.clustatis$param$Noise_cluster
    gpmax=res.clustatis$param$gpmax
    if (Noise_cluster==FALSE)
    {
      dif=res.clustatis$diff_crit_ng[1,]
      quality=res.clustatis$overall_homogeneity_ng[1,]
    }else{
      dif=res.clustatis$diff_crit_ng
      quality=res.clustatis$overall_homogeneity_ng
    }

    dev.new()
    barplot(dif,xlab="Nb clusters",ylab="delta", main="Variation of criterion before consolidation",
            axisnames = TRUE,names.arg = paste((gpmax):2,"->",(gpmax-1):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

    dev.new()
    barplot(quality,xlab="Nb clusters",ylab="Overall homogeneity (%)", main="Overall homogeneity before consolidation (%)",
            axisnames = TRUE,names.arg = paste((gpmax):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

    if (Noise_cluster==FALSE)
    {
      dif=res.clustatis$diff_crit_ng[2,]
      quality=res.clustatis$overall_homogeneity_ng[2,]

      dev.new()
      barplot(dif,xlab="Nb clusters",ylab="delta", main="Variation of criterion after consolidation",
              axisnames = TRUE,names.arg = paste((gpmax):2,"->",(gpmax-1):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")

      dev.new()
      barplot(quality,xlab="Nb clusters",ylab="Overall homogeneity (%)", main="Overall homogeneity after consolidation (%)",
              axisnames = TRUE,names.arg = paste((gpmax):1),las=2,cex.names = 0.6,cex.main=1.2,col = "blue")
    }
  }

  if(Graph_groups==TRUE)
  {
    for (i in 1:ngroups)
    {
      if(res.clustatis$type=="H+C")
      {
        C=res.clustatis[[ngroups]]$coord[[i]]
        pouriner=res.clustatis[[ngroups]]$inertia[[i]]
      }else{
        C=res.clustatis$coord[[i]]
        pouriner=res.clustatis$inertia[[i]]
      }

      un=axes[1]
      deux=axes[2]
      dev.new()
      plot(C[,un],C[,deux],type="n",lwd=5,pch=16,xlab=paste("Dim",axes[1], "(" ,pouriner[un],"%)"), ylab=paste("Dim",axes[2], "(" ,pouriner[deux],"%)"),xlim=c(min(C[,un])-0.2,max(C[,un])+0.2), ylim=c(min(C[,deux])-0.1,max(C[,deux])+0.1))
      text(C[,un],C[,deux],rownames(C),col=col, cex=cex, font = font)
      abline(h=0,v=0)
      title(paste("Cluster",i))
    }
  }

  if(Graph_weights==TRUE)
  {
    if(res.clustatis$type=="H+C")
    {
      res.clustatis=res.clustatis[[ngroups]]
    }
    for (i in 1:ngroups)
    {
      dev.new()
      barplot(res.clustatis$weights[[i]])
      title(paste("Weights in the cluster",i))
    }
  }

}
