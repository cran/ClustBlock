##=============================================================================
##' @title Print the ClusMB or clustering on STATIS axes results
##'
##' @usage
##' \method{print}{clusRows}(x, ...)
##'
##' @description
##' Print the ClusMB or clustering on STATIS axes results
##'
##' @param x object of class 'clusRows'
##'
##' @param ... further arguments passed to or from other methods
##'
##' @keywords quantitative CATA RATA
##'
##' @seealso   \code{\link{ClusMB}}, \code{\link{clustRowsOnStatisAxes}}
##'
##' @export


##=============================================================================

print.clusRows=function(x, ...)
{
  res.ClusRows=x
  if(inherits(res.ClusRows, "clusRows")==FALSE)
  {
    stop("The class of the object must be 'clusRows'")
  }

  cat(paste("Clustering of rows in Multi-Block context","\n" ))
  cat(paste("number of blocks:",res.ClusRows$param$nblo, "\n"))
  cat(paste("number of rows:",res.ClusRows$param$n, "\n"))
  cat(paste("Method used:",res.ClusRows$type, "\n"))
}
