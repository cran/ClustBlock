## =============================================================================

##' @title Displays the ClusMB and clustRowsOnstatisAxes graphs
##'
##' @usage
##'  \method{plot}{clusRows}(x, ...)
##'
##' @description
##' This function plots the dendrogram of ClusMB or clustRowsOnstatisAxes
##'
##'
##' @param x object of class 'clusRows'
##'
##' @param ... further arguments passed to or from other methods
##'
##
##' @return the dendrogram
##'
##'
##' @keywords quantitative CATA RATA
##'
##' @importFrom graphics legend
##'
##' @examples
##' ##'
##' #####projective mapping####
##' library(ClustBlock)
##' data(smoo)
##' res1=ClusMB(smoo, rep(2,24))
##' plot(res1)
##'
##' @seealso   \code{\link{ClusMB}}, \code{\link{clustRowsOnStatisAxes}}
##'
##' @export
##'

## =============================================================================




plot.clusRows <- function(x, ...) {
  Sys.setenv("_R_CHECK_SYSTEM_CLOCK_" = 0)

  res.clusRows <- x
  if (inherits(res.clusRows, "clusRows") == FALSE) {
    stop("The class of the object must be 'clusRows'")
  }


  # show graphical representation
  dev.new()
  plot(res.clusRows$dend, hang = -1, main = "Dendrogram", axes = TRUE, ylab = "Height", sub = "", xlab = "", ...)
}
