
##=============================================================================

##' @title Print the CLUSCATA results
##'
##'
##' @usage
##' \method{print}{cluscata}(x, ...)
##'
##' @description
##' Print the CLUSCATA results
##'
##'
##' @param x object of class 'cluscata'
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @keywords CATA
##'
##' @seealso   \code{\link{cluscata}} , \code{\link{cluscata_kmeans}}
##'
##' @export


##=============================================================================


print.cluscata=function(x, ...)
{
  res.cluscata=x
  if(inherits(res.cluscata, "cluscata")==FALSE)
  {
    stop("The class of the object must be 'cluscata'")
  }

  if (res.cluscata$type=="H+C")
  {
    cat("Hierarchical clustering of binary blocks with consolidation \n")
    cat(paste("number of subjects:",res.cluscata$param$nblo, "\n"))
    cat(paste("number of products by subject:",res.cluscata$param$n, "\n"))
    cat(paste("number of attributes by subject:",res.cluscata$param$nvar, "\n"))
    cat(paste0("consolidation for K in 1:", res.cluscata$param$gpmax, "\n"))
    cat(paste("noise Cluster:", res.cluscata$param$Noise_cluster, "\n"))
    cat("\n" )
    cat("\n" )
    cat("$partitionK or [[K]]: results with K clusters after consolidation \n")
    cat("$cutree_k$partitionK: partition in K clusters before consolidation \n")
  }else if (res.cluscata$type=="K"){
    cat("Partitionning clustering of binary Blocks\n")
    cat(paste("number of subjects:", res.cluscata$param$nblo, "\n"))
    cat(paste("number of products by subject:",res.cluscata$param$n, "\n"))
    cat(paste("number of attributes by subject:",res.cluscata$param$nvar, "\n"))
    cat(paste("number of clusters:" ,res.cluscata$param$ng, "\n"))
    cat(paste("threshold for Noise Cluster:", res.cluscata$rho, "\n"))
  }
}
