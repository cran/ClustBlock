##=============================================================================


##' @title Print the CATATIS results
##'
##' @usage
##' \method{print}{catatis}(x, ...)
##'
##' @description
##' Print the CATATIS results
##'
##'
##' @param x object of class 'catatis'
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##'
##' @keywords CATA
##'
##' @seealso   \code{\link{catatis}}
##'
##' @export


##=============================================================================




print.catatis=function(x, ...)
{
  res.catatis=x
  if(class(res.catatis)!="catatis")
  {
    stop("The class of the object must be 'catatis'")
  }

  cat("CATATIS method on binary blocks\n")
  cat(paste("number of subjects:",res.catatis$param$nblo, "\n"))
  cat(paste("number of products by subject:",res.catatis$param$n, "\n"))
  cat(paste("number of attributes by subject:",res.catatis$param$nvar, "\n"))
}
