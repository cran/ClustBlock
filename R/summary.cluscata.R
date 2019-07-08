##=============================================================================


##' @title Show the CLUSCATA results
##'
##' @usage
##' \method{summary}{cluscata}(object, ngroups=NULL, ...)
##'
##' @description
##' This function shows the cluscata results
##'
##'
##' @param object object of class 'cluscata'.
##'
##' @param ngroups number of groups to consider. Ignored for cluscata_kmeans results. Default: recommended  number of clusters
##'
##' @param ... further arguments passed to or from other methods
##'
##'
##' @return the CLUSCATA principal results
##'
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item weights: weight associated with each subject in its cluster
##'          \item rho: the threshold for the noise cluster
##'          \item test_one_cluster: decision and pvalue to know if there is more than one cluster
##'          }
##'
##' @keywords CATA
##'
##' @seealso   \code{\link{cluscata}} , \code{\link{cluscata_kmeans}}
##'
##' @export


##=============================================================================





summary.cluscata=function(object, ngroups=NULL, ...)
{
  res.cluscata=object
  if(class(res.cluscata)!="cluscata")
  {
    stop("The class of the object must be 'cluscata'")
  }

  if(is.null(ngroups) | res.cluscata$type=="K")
  {
    ngroups=res.cluscata$param$ng
  }

  if(res.cluscata$type=="H+C")
  {
    if(ngroups>res.cluscata$param$gpmax)
    {
      stop("ngroups>gpmax")
    }
  }

  if(res.cluscata$type=="H+C")
  {
    test_one_cluster=res.cluscata$test_one_cluster
    res.cluscata=res.cluscata[[ngroups]]
  }else{
    test_one_cluster="No test"
  }
  NameBlocks=rownames(res.cluscata$group)

  liste_groups=NULL
  for (i in 1:ngroups)
  {
    liste_groups[[i]]=NameBlocks[res.cluscata$group==i]
  }
  names(liste_groups)=paste("Cluster", 1:ngroups)
  if (res.cluscata$rho>0)
  {
    liste_groups[[ngroups+1]]=NameBlocks[res.cluscata$group=="K+1"]
    names(liste_groups)[ngroups+1]="Noise cluster (K+1)"
  }

  res=list(groups=liste_groups, homogeneity=res.cluscata$homogeneity,
           weights=res.cluscata$weights, rho=res.cluscata$rho, test_one_cluster=test_one_cluster)
  return(res)
}
