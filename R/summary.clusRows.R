##=============================================================================
##' @title Show the ClusMB or clustering on STATIS axes results
##'
##' @usage
##' \method{summary}{clusRows}(object, ...)
##'
##' @description
##' This function shows the ClusMB or clustering on STATIS axes results
##'
##'
##' @param object object of class 'clusRows'.
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item groups: clustering partition
##'          \item nbClustRetained: the number of clusters retained
##'          \item nbgH: Advised number of clusters per Hartigan index
##'          \item nbgCH: Advised number of clusters per Calinski-Harabasz index
##'          }
##'
##'
##' @keywords quantitative CATA RATA
##'
##'
##' @seealso   \code{\link{ClusMB}}, \code{\link{clustRowsOnStatisAxes}}
##'
##' @export


##=============================================================================


summary.clusRows=function(object, ...)
{
  res.ClusRows=object
  if(inherits(res.ClusRows, "clusRows")==FALSE)
  {
    stop("The class of the object must be 'clusRows'")
  }


  res=list(groups=res.ClusRows$group, nbClustRetained= max(res.ClusRows$group),
           nbgH=res.ClusRows$nbgH, nbgCH=res.ClusRows$nbgCH)
  return(res)
}
