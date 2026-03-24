## =============================================================================
##' @title Print the CLUSCATA-liking results
##'
##' @usage
##' \method{print}{cluscata_liking}(x, ...)
##'
##' @description
##'  Print the CLUSCATA-liking results
##'
##' @param x object of class 'cluscata_liking'
##'
##' @param ... further arguments passed to or from other methods
##'
##' @keywords CATA liking
##'
##' @seealso   \code{\link{cluscata_liking}}
##'
##' @export


## =============================================================================

print.cluscata_liking <- function(x, ...) {
  res.cluscata_liking <- x
  if (inherits(res.cluscata_liking, "cluscata_liking") == FALSE) {
    stop("The class of the object must be 'cluscata_liking'")
  }

  cat(paste("CLUSCATA-liking", "\n"))
  cat(paste("number of subjects:", res.cluscata_liking$param$nblo, "\n"))
  cat(paste("number of products:", res.cluscata_liking$param$n, "\n"))
  cat(paste("number of attributes:", res.cluscata_liking$param$nvar, "\n"))
}
