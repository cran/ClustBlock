##=============================================================================


##' @title Displays the STATIS graphs
##'
##'
##' @usage
##' \method{plot}{statis}(x, axes=c(1,2), Graph_obj=TRUE,
##' Graph_weights=TRUE, tit="STATIS", col=NULL, cex=1, ...)
##'
##' @description
##' This function plots the STATIS map and STATIS weights
##'
##'
##' @param x object of class 'statis'
##'
##' @param Graph_obj logical. Should the compromise graphical representation be plotted? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##' @param axes  numerical vector (length 2). Axes to be ploted. Default: c(1,2)
##'
##' @param tit string. Title for the objects graphical representation. Default: 'STATIS'
##'
##' @param col vector. Color for each object. If NULL, col=rainbow(nrow(Data)). Default: NULL
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return the STATIS graphs
##'
##' @keywords quantitative
##'
##'
##' @examples
##'
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  st=statis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks)
##'  plot(st, axes=c(1,3), Graph_weights=FALSE)
##'
##' @seealso \code{\link{statis}}
##'
##' @export

##=============================================================================



plot.statis=function(x, axes=c(1,2), Graph_obj=TRUE, Graph_weights=TRUE, tit="STATIS", col=NULL, cex=1, ...)
{

  res.statis=x
  if(class(res.statis)!="statis")
  {
    stop("The class of the object must be 'statis'")
  }

  if(Graph_obj==TRUE)
  {
    #show graphical representation
    C=res.statis$coord
    if(is.null(col)==TRUE)
    {
      col=rainbow(nrow(C))
    }
    pouriner=res.statis$inertia
    un=axes[1]
    deux=axes[2]
    dev.new()
    plot(C[,un],C[,deux],type="n",lwd=5,pch=16,xlab=paste("Dim",axes[1], "(" ,pouriner[un],"%)"), ylab=paste("Dim",axes[2], "(" ,pouriner[deux],"%)"),xlim=c(min(C[,un])-0.2,max(C[,un])+0.2), ylim=c(min(C[,deux])-0.2,max(C[,deux])+0.2))
    text(C[,un],C[,deux],rownames(C),col=col, cex=cex)
    abline(h=0,v=0)
    title(tit)


    #projection of each object of each block

    dev.new()
    n=res.statis$param$n
    nblo=res.statis$param$nblo
    configs=res.statis$proj_config
    plot(C[,un],C[,deux],type="n",lwd=5,pch=16,xlab=paste("Dim",axes[1], "(" ,pouriner[un],"%)"), ylab=paste("Dim",axes[2], "(" ,pouriner[deux],"%)"),xlim=c(min(C[,un])-0.25,max(C[,un])+0.25), ylim=c(min(C[,deux])-0.25,max(C[,deux])+0.25))
    text(C[,un],C[,deux],rownames(C),col=rainbow(n),font=2)
    for (l in 1:nblo)
    {
      points(configs[,un,l],configs[,deux,l],col=rainbow(n),pch=20)
      arrows(C[,un],C[,deux],configs[,un,l],configs[,deux,l],col=rainbow(n),lwd=0.5,lty = 2,angle=0)
    }
    abline(h=0,v=0)
    title(tit)
  }

  if (Graph_weights==TRUE)
  {
    dev.new()
    barplot(res.statis$weights)
    title(paste("Weights"))
  }


}
