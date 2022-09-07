##=============================================================================


##' @title Show the CATATIS results
##'
##' @usage
##' \method{summary}{catatis}(object, ...)
##'
##' @description
##' This function shows the CATATIS results
##'
##'
##' @param object object of class 'catatis'.
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item homogeneity: homogeneity of the subjects (in percentage)
##'          \item weights: the weights associated with the subjects to build the compromise
##'          \item eigenvalues: the eigenvalues associated to the correspondance analysis
##'          \item inertia: the percentage of total variance explained by each axis of the CA
##'          }
##'
##'
##' @keywords CATA
##'
##'
##' @seealso   \code{\link{catatis}}
##'
##' @export


##=============================================================================





summary.catatis=function(object, ...)
{
  res.catatis=object
  if(inherits(res.catatis, "catatis")==FALSE)
  {
    stop("The class of the object must be 'catatis'")
  }

  NameBlocks=names(res.catatis$weights)


  res=list(homogeneity=res.catatis$homog, weights=res.catatis$weights,
           eigenvalues=res.catatis$eigenvalues, inertia=res.catatis$inertia)
  return(res)
}
