## =============================================================================
##' @title Combination of CATA and liking data for CLUSCATA-liking
##'
##' @description
##' For CLUSCATA-liking, this preprocessing is needed.
##'
##' @usage
##' combinCATALiking(cata, liking, center=TRUE, scale=FALSE)
##'
##'
##' @param cata data frame or matrix where the blocks of binary variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param liking data frame or matrix where the products are in rows and the assessors in columns
##'
##' @param center Centering of consumer liking. Default: TRUE
##'
##' @param scale Scaling  of consumer liking. Default: FALSE
##'
##' @return Combined data
##'
##' @keywords CATA liking
##'
##' @references
##' Vigneau, E., Cariou, V., Giacalone, D., Berget, I., & Llobell, F. (2022). Combining hedonic information and CATA description for consumer segmentation. Food Quality and Preference, 95, 104358.
##'
##'
##' @examples
##' data(cata_ryebread)
##' data(liking_ryebread)
##' cataliking=combinCATALiking(cata_ryebread, liking_ryebread)
##'
##' @seealso   \code{\link{cluscata_liking}}
##'
##' @export
## =============================================================================


combinCATALiking=function(cata, liking, center=TRUE, scale=FALSE)
{
  cata=as.matrix(cata)
  likingc=scale(liking, center=center, scale=scale)
  p=ncol(liking)
  Q=ncol(cata)/p
  n=nrow(liking)
  Blocks=rep(Q, p)
  J = rep(1:p, times = Blocks)
  catalM=matrix(NA,n,p*Q)
  for (j in 1:p) {
    yj=likingc[,j]
    Xj=cata[, j==J]
    catalM[,j==J]=diag(yj)%*%Xj
  }
  colnames(catalM)=colnames(cata)
  rownames(catalM)=rownames(cata)
  return(catalM)
}
