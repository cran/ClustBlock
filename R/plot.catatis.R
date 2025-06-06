## =============================================================================

##' @title Displays the CATATIS graphs
##'
##' @usage
##'  \method{plot}{catatis}(x, Graph=TRUE, Graph_weights=TRUE, Graph_eig=TRUE,
##'   axes=c(1,2), tit="CATATIS", cex=1, col.obj="blue", col.attr="red", ...)
##'
##' @description
##' This function plots the CATATIS map and CATATIS weights
##'
##'
##' @param x object of class 'catatis'
##'
##' @param Graph logical. Show the graphical representation? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##' @param Graph_eig logical. Should the barplot of the eigenvalues be plotted? Only with Graph=TRUE. Default: TRUE
##'
##' @param axes  numerical vector (length 2). Axes to be plotted
##'
##' @param tit string. Title for the graphical representation. Default: 'CATATIS'
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param col.obj numerical or string. Color for the objects points. Default: "blue"
##'
##' @param col.attr numerical or string. Color for the attributes points. Default: "red"
##'
##' @param ... further arguments passed to or from other methods
##'
##
##' @importFrom FactoMineR plot.CA
##'
##' @return the CATATIS map
##'
##'
##' @keywords CATA RATA
##'
##' @examples
##'  \donttest{
##' data(straw)
##' res.cat=catatis(straw, nblo=114)
##' plot(res.cat, Graph_weights=FALSE, axes=c(1,3))
##' }
##'
##' @seealso   \code{\link{catatis}}
##'
##' @export
##'

## =============================================================================




plot.catatis <- function(x, Graph = TRUE, Graph_weights = TRUE, Graph_eig = TRUE, axes = c(1, 2), tit = "CATATIS", cex = 1,
                         col.obj = "blue", col.attr = "red", ...) {
  res.catatis <- x
  if (inherits(res.catatis, "catatis") == FALSE) {
    stop("The class of the object must be 'catatis'")
  }


  # show graphical representation
  if (Graph == TRUE) {
    vp <- res.catatis$CA$eig[, 1]
    if (Graph_eig == TRUE) {
      dev.new()
      barplot(vp, col = "blue", main = "Eigenvalues")
    }
    dev.new()
    options(ggrepel.max.overlaps = Inf)
    print(plot.CA(res.catatis$CA, axes = axes, title = tit, cex = cex, col.row = col.obj, col.col = col.attr))
  }

  if (Graph_weights == TRUE) {
    dev.new()
    barplot(res.catatis$weights)
    title(paste("Weights"))
  }
}
