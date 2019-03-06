##=============================================================================


##' @title Show the CLUSTATIS results
##'
##' @usage
##' \method{summary}{clustatis}(object, ngroups=NULL, ...)
##'
##' @description
##' This function shows the clustatis results
##'
##'
##' @param object object of class 'clustatis'.
##'
##' @param ngroups number of groups to consider. Ignored for clustatis_kmeans results. Default: recommended  number of clusters
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return the CLUSTATIS principal results
##'
##'
##'
##' @keywords quantitative
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item weights: weight associated with each block in its cluster
##'          \item rho: the threshold for the noise cluster
##'          }
##'
##'
##' @seealso   \code{\link{clustatis}} , \code{\link{clustatis_kmeans}}
##'
##' @export


##=============================================================================





summary.clustatis=function(object, ngroups=NULL, ...)
{
  res.clustatis=object
  if(class(res.clustatis)!="clustatis")
  {
    stop("The class of the object must be 'clustatis'")
  }

  if(is.null(ngroups) | res.clustatis$type=="K")
  {
    ngroups=res.clustatis$param$ng
  }

  if(res.clustatis$type=="H+C")
  {
    if(ngroups>res.clustatis$param$gpmax)
    {
      stop("ngroups>gpmax")
    }
  }

  if(res.clustatis$type=="H+C")
  {
    res.clustatis=res.clustatis[[ngroups]]
  }
  NameBlocks=rownames(res.clustatis$group)

  liste_groups=NULL
  for (i in 1:ngroups)
  {
    liste_groups[[i]]=NameBlocks[res.clustatis$group==i]
  }
  names(liste_groups)=paste("Cluster", 1:ngroups)
  if (res.clustatis$rho>0)
  {
    liste_groups[[ngroups+1]]=NameBlocks[res.clustatis$group=="K+1"]
    names(liste_groups)[ngroups+1]="Noise cluster (K+1)"
  }

  res=list(groups=liste_groups, homogeneity=res.clustatis$homogeneity,
        weights=res.clustatis$weights, rho=res.clustatis$rho)
  return(res)
}
