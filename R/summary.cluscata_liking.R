## =============================================================================
##' @title Show the CLUSCATA-liking results
##'
##' @usage
##' \method{summary}{cluscata_liking}(object, ngroups=NULL, ...)
##'
##' @description
##' This function shows the cluscata_liking results
##'
##'
##' @param object object of class 'cluscata_liking'.
##'
##' @param ngroups number of groups to consider. Default: recommended  number of clusters
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return the CLUSCATA-liking principal results
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition
##'          \item groupsvector: the groups in a vector
##'          }
##'
##' @keywords CATA liking
##'
##' @seealso   \code{\link{cluscata_liking}} , \code{\link{plot.cluscata_liking}}
##'
##' @export


## =============================================================================





summary.cluscata_liking <- function(object, ngroups = NULL, ...) {
  res.cluscata <- object
  if (inherits(res.cluscata, "cluscata_liking") == FALSE) {
    stop("The class of the object must be 'cluscata_liking'")
  }

  if (is.null(ngroups) | res.cluscata$type == "K") {
    ngroups <- res.cluscata$param$ng
  }

  if (res.cluscata$type == "H+C") {
    if (ngroups > res.cluscata$param$gpmax) {
      stop("ngroups>gpmax")
    }
    res.cluscata <- res.cluscata[[ngroups]]
  }


  NameBlocks <- rownames(res.cluscata$group)

  liste_groups <- list()
  for (i in 1:ngroups)
  {
    liste_groups[[i]] <- NameBlocks[res.cluscata$group == i]
  }
  names(liste_groups) <- paste("Cluster", 1:ngroups)

  res <- list(
    groups = liste_groups, groupsvector = res.cluscata$group)
  return(res)
}
