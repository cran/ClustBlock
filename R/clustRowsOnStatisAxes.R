##' @title Perform a cluster analysis of rows in a Multi-block context with clustering on STATIS axes
##'
##' @description
##' Clustering of rows (products in sensory analysis) in a Multi-block context.
##' The STATIS method is followed by a hierarchical algorithm.
##'
##' @usage
##' clustRowsOnStatisAxes(Data, Blocks, NameBlocks=NULL, scale=FALSE,
##' nclust=NULL, gpmax=6, ncomp=5)
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data.
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##' @param nclust numerical. Number of clusters to consider. If NULL, the Hartigan index advice is taken.
##'
##' @param gpmax logical. What is maximum number of clusters to consider?  min(6, number of blocks -2)
##'
##' @param ncomp numerical. Number of axes to consider. Default:5
##'
##' @return
##'         \itemize{
##'          \item group: the clustering partition.
##'          \item nbgH: Advised number of clusters per Hartigan index
##'          \item nbgCH: Advised number of clusters per Calinski-Harabasz index
##'          \item cutree_k: the partition obtained by cutting the dendrogram in K clusters
##'          \item dend: The dendrogram
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##' @keywords quantitative CATA RATA
##'
##' @references
##' Llobell, F., & Giacalone, D. (2025). Two Methods for Clustering Products in a Sensory Study: STATIS and ClusMB. Journal of Sensory Studies, 40(1), e70024.\cr
##'
##' @importFrom  stats dist hclust kmeans cutree
##'
##'
##' @examples
##'
##' #####projective mapping####
##' library(ClustBlock)
##' data(smoo)
##' res1=clustRowsOnStatisAxes(smoo, rep(2,24))
##' summary(res1)
##' indicesClusters(smoo, rep(2,24), res1$group)
##'
##' ####CATA####
##' data(fish)
##' Data=fish[1:66,2:30]
##' chang2=change_cata_format2(Data, nprod= 6, nattr= 27, nsub = 11, nsess= 1)
##' res2=clustRowsOnStatisAxes(Data= chang2$Datafinal, Blocks= rep(27, 11))
##' indicesClusters(Data= chang2$Datafinal, Blocks= rep(27, 11),cut = res2$group, center=FALSE)
##'
##' @seealso   \code{\link{indicesClusters}}, \code{\link{summary.clusRows}} , \code{\link{ClusMB}}
##'
##' @export


## =============================================================================


clustRowsOnStatisAxes <- function(Data, Blocks, NameBlocks = NULL, scale = FALSE,
                                  nclust = NULL, gpmax = 6, ncomp = 5) {
  nblo <- length(Blocks)
  n <- nrow(Data)
  if (is.null(NameBlocks)) NameBlocks <- paste("B", 1:nblo, sep = "-")
  if (is.null(rownames(Data))) rownames(Data) <- paste0("X", 1:nrow(Data))
  if (is.null(colnames(Data))) colnames(Data) <- paste0("Y", 1:ncol(Data))

  # parapet for numerical Data
  for (i in 1:ncol(Data))
  {
    if (is.numeric(Data[, i]) == FALSE) {
      stop(paste("The data must be numeric (column", i, ")"))
    }
  }

  # parapet for number of objects
  if (n < 3) {
    stop("At least 3 objects are required")
  }

  # parapet for number of blocks
  if (nblo < 2) {
    stop("At least 2 blocks are required")
  }


  # parapet for Blocks
  if (sum(Blocks) != ncol(Data)) {
    stop("Error with Blocks")
  }

  # Parapet for NameBlocks
  if (length(NameBlocks) != nblo) {
    stop("Error with the length of NameBlocks")
  }


  # parapet for scale: no constant variable
  if (scale == TRUE) {
    for (i in 1:ncol(Data))
    {
      if (sd(Data[, i]) == 0) {
        stop(paste("Column", i, "is constant"))
      }
    }
  }

  # no NA
  if (sum(is.na(Data)) > 0) {
    print("NA detected:")
    tabna <- which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }


  st <- statis(Data, Blocks, Graph_obj = FALSE, Graph_weights = FALSE, scale = scale)

  res2 <- hclust(dist(st$coord[, 1:ncomp]), method = "ward.D2")
  res2$height <- (res2$height / sqrt(2))**2
  dev.new()
  plot(res2, hang = -1, main = "Dendrogram", axes = TRUE, ylab = "Height", sub = "", xlab = "")

  # number of clusters advised by hartigan
  criter <- sort(cumsum(res2$height), decreasing = TRUE)
  H <- NULL
  for (k in 1:min(gpmax, n - 2))
  {
    H[k] <- (criter[k] / criter[k + 1] - 1) * (n - k - 1)
  }
  nbgroup_hart <- which.max(H[-(min(gpmax, n - 2))] - H[-1]) + 1
  cat(paste("Recommended number of clusters =", nbgroup_hart), "\n")

  # number of clusters advised by Calinsi Harabasz
  bgss <- NULL
  for (i in 1:gpmax)
  {
    bgss[i] <- nblo - criter[i]
  }

  CH <- NULL
  CH[1] <- 0
  for (k in 2:gpmax)
  {
    num <- (bgss[k] / (k - 1))
    den <- (criter[k]) / (n - k)
    CH[k] <- num / den
  }
  nbgroup_ch <- which.max(CH)

  # take the number of clusters suggested by Hartigan if not privided
  if (is.null(nclust)) {
    nclust <- nbgroup_hart
  }

  clust <- cutree(res2, nclust)

  cutree_k <- list()
  for (K in 1:gpmax)
  {
    coupe <- cutree(res2, K)
    cutree_k[[K]] <- coupe
  }
  names(cutree_k) <- paste0("partition", 1:gpmax)



  resfin <- list(
    group = clust, nbgH = nbgroup_hart,
    nbgCH = nbgroup_ch, cutree_k = cutree_k, dend = res2,
    param = list(nblo = nblo, gpmax = gpmax, n = n),
    type = "clusonstatisaxes"
  )

  class(resfin) <- "clusRows"

  return(resfin)
}
