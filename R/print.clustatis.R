## =============================================================================

##' @title Print the CLUSTATIS results
##' @usage
##' \method{print}{clustatis}(x, ...)
##'
##' @description
##' Print the CLUSTATIS results
##'
##'
##' @param x object of class 'clustatis'
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##'
##' @keywords quantitative
##'
##' @seealso   \code{\link{clustatis}}  , \code{\link{clustatis_kmeans}}
##'
##' @export


## =============================================================================


print.clustatis <- function(x, ...) {
  res.clustatis <- x
  if (inherits(res.clustatis, "clustatis") == FALSE) {
    stop("The class of the object must be 'clustatis'")
  }

  if (res.clustatis$type == "H+C") {
    cat("Hierarchical clustering of quantitative blocks with consolidation \n")
    cat(paste("number of blocks:", res.clustatis$param$nblo, "\n"))
    cat(paste("number of objects:", res.clustatis$param$n, "\n"))
    cat(paste0("consolidation for K in 1:", res.clustatis$param$gpmax, "\n"))
    cat(paste("noise Cluster:", res.clustatis$param$Noise_cluster, "\n"))
    cat("\n")
    cat("\n")
    cat("$partitionK or [[K]]: results with K clusters after consolidation \n")
    cat("$cutree_k$partitionK: partition in K clusters before consolidation \n")
  } else if (res.clustatis$type == "K") {
    cat("partitioning clustering of quantitative Blocks\n")
    cat(paste("number of blocks:", res.clustatis$param$nblo, "\n"))
    cat(paste("number of objects:", res.clustatis$param$n, "\n"))
    cat(paste("number of clusters:", res.clustatis$param$ng, "\n"))
    for (i in 1:res.clustatis$param$ng)
    {
      cat(paste("threshold for Noise Cluster", i, ":", res.clustatis$rho[i], "\n"))
    }
  }
}
