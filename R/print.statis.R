##=============================================================================


##' @title Print the STATIS results
##'
##' @usage
##' \method{print}{statis}(x, ...)
##'
##' @description
##' Print the STATIS results
##'
##'
##' @param x object of class 'statis'
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##'
##' @keywords quantitative
##'
##' @seealso   \code{\link{statis}}
##'
##' @export


##=============================================================================




print.statis=function(x, ...)
{
  res.statis=x
  if(inherits(res.statis, "statis")==FALSE)
  {
    stop("The class of the object must be 'statis'")
  }

    cat("STATIS method on quantitative blocks\n")
    cat(paste("number of blocks:",res.statis$param$nblo, "\n"))
    cat(paste("number of objects:",res.statis$param$n, "\n"))
}
