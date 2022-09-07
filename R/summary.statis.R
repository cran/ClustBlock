##=============================================================================


##' @title Show the STATIS results
##'
##' @usage
##' \method{summary}{statis}(object, ...)
##'
##' @description
##' This function shows the STATIS results
##'
##'
##' @param object object of class 'statis'.
##'
##' @param ... further arguments passed to or from other methods
##'
##
##'
##'
###' @return a list with:
##'         \itemize{
##'          \item homogeneity: homogeneity of the blocks (in percentage)
##'          \item weights: the weights associated with the blocks to build the compromise
##'          \item eigenvalues: the eigenvalues of the svd decomposition
##'          \item inertia: the percentage of total variance explained by each axis
##'          }
##'
##'
##' @keywords quantitative
##'
##' @seealso   \code{\link{statis}}
##'
##' @export


##=============================================================================





summary.statis=function(object, ...)
{
  res.statis=object
  if(inherits(res.statis, "statis")==FALSE)
  {
    stop("The class of the object must be 'statis'")
  }

  NameBlocks=names(res.statis$weights)


  res=list(homogeneity=res.statis$homogeneity, weights=res.statis$weights,
            eigenvalues=res.statis$eigenvalues, inertia=res.statis$inertia)
  return(res)
}
