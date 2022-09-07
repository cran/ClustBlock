##=============================================================================


##' @title Displays the STATIS graphs
##'
##'
##' @usage
##' \method{plot}{statis}(x, axes=c(1,2), Graph_obj=TRUE,
##' Graph_weights=TRUE, Graph_eig=TRUE, tit="STATIS", col=NULL, cex=1, font=1,
##' xlim=NULL, ylim=NULL, ...)
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
##' @param Graph_eig logical. Should the barplot of the eigenvalues be plotted? Only with Graph_obj=TRUE. Default: TRUE
##'
##' @param axes  numerical vector (length 2). Axes to be plotted. Default: c(1,2)
##'
##' @param tit string. Title for the objects graphical representation. Default: 'STATIS'
##'
##' @param col vector. Color for each object. If NULL, col=rainbow(nrow(Data)). Default: NULL
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param font numerical. Integer specifying font to use for text. 1=plain, 2=bold, 3=italic, 4=bold italic, 5=symbol. Default: 1
##'
##' @param xlim numerical vector (length 2). Minimum and maximum for x coordinates.
##'
##' @param ylim numerical vector (length 2). Minimum and maximum for y coordinates.
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
##' \donttest{
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  st=statis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks)
##'  plot(st, axes=c(1,3), Graph_weights=FALSE)
##'
##'  }
##'
##' @seealso \code{\link{statis}}
##'
##' @export

##=============================================================================



plot.statis=function(x, axes=c(1,2), Graph_obj=TRUE, Graph_weights=TRUE, Graph_eig=TRUE,
                     tit="STATIS", col=NULL, cex=1, font=1, xlim=NULL, ylim=NULL, ...)
{

  res.statis=x
  if(inherits(res.statis, "statis")==FALSE)
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
    vp=res.statis$eigenvalues
    if (Graph_eig==TRUE)
    {
      dev.new()
      barplot(vp, col="blue", main="Eigenvalues")
    }
    pouriner=res.statis$inertia
    un=axes[1]
    deux=axes[2]
    if(is.null(xlim)==TRUE)
    {
      xlim=c(min(C[,un])-0.2,max(C[,un])+0.2)
    }
    if(is.null(ylim)==TRUE)
    {
      ylim=c(min(C[,deux])-0.2,max(C[,deux])+0.2)
    }
    dev.new()
    plot(C[,un],C[,deux],type="n",lwd=5,pch=16,xlab=paste("Dim",axes[1], "(" ,pouriner[un],"%)"), ylab=paste("Dim",axes[2], "(" ,pouriner[deux],"%)"), xlim=xlim, ylim=ylim)
    text(C[,un],C[,deux],rownames(C),col=col, cex=cex, font=font)
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
