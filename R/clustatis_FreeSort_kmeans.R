## =============================================================================


##' @title Compute the CLUSTATIS partitioning algorithm on free sorting data
##'
##' @description
##' partitioning algorithm for Free Sorting data. Each cluster is associated with a compromise
##' computed by the STATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##' clustatis_FreeSort_kmeans(Data, NameSub=NULL, clust, nstart=100, rho=0,Itermax=30,
##' Graph_groups=TRUE, Graph_weights=FALSE,  print_attempt=FALSE)
##'
##'
##' @param Data data frame or matrix. Corresponds to all variables that contain subjects results. Each column corresponds to a subject and gives the groups to which the products (rows) are assigned
##'
##' @param NameSub string vector. Name of each subject. Length must be equal to the number of clumn of the Data. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param clust numerical vector or integer. Initial partition or number of starting partitions if integer. If numerical vector, the numbers must be 1,2,3,...,number of clusters
##'
##' @param nstart integer. Number of starting partitions. Default: 100
##'
##' @param rho numerical or vector between 0 and 1. Threshold for the noise cluster. Default:0. If you want a different threshold for each cluster, you can provide a vector.
##'
##' @param Itermax numerical. Maximum of iterations by partitioning algorithm. Default: 30
##'
##' @param Graph_groups logical. Should each cluster compromise be plotted? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights in each cluster be plotted? Default: FALSE
##'
##' @param print_attempt logical. Print the number of remaining attempts in the multi-start case? Default: FALSE
##'
##'
##'
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition. If rho>0, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold(s) for the noise cluster
##'          \item homogeneity: percentage of homogeneity of the subjects in each cluster and the overall homogeneity
##'          \item rv_with_compromise: RV coefficient of each subject with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item comp_RV: RV coefficient between the compromises associated with the various clusters
##'          \item compromise: the W compromise of each cluster
##'          \item coord: the coordinates of objects of each cluster
##'          \item inertia: percentage of total variance explained by each axis for each cluster
##'          \item rv_all_cluster: the RV coefficient between each subject and each cluster compromise
##'          \item criterion: the CLUSTATIS criterion error
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##'
##' @keywords FreeSorting
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods. Application to sensometrics. Food Quality and Preference, in Press.\cr
##' Llobell, F., Vigneau, E., Qannari, E. M. (2019). Clustering datasets by means of CLUSTATIS with identification of atypical datasets. Application to sensometrics. Food Quality and Preference, 75, 97-104.
##'
##'
##'
##' @examples
##' data(choc)
##' res.clu=clustatis_FreeSort_kmeans(choc, clust=2)
##' plot(res.clu, Graph_groups=FALSE, Graph_weights=TRUE)
##' summary(res.clu)
##'
##' @seealso   \code{\link{clustatis_FreeSort}}, \code{\link{preprocess_FreeSort}}, \code{\link{summary.clustatis}}, , \code{\link{plot.clustatis}}
##'
##' @export


## =============================================================================


clustatis_FreeSort_kmeans <- function(Data, NameSub = NULL, clust, nstart = 100, rho = 0,
                                      Itermax = 30, Graph_groups = TRUE,
                                      Graph_weights = FALSE, print_attempt = FALSE) {
  prepro <- preprocess_FreeSort(Data, NameSub = NameSub)

  a <- clustatis_kmeans(
    Data = prepro$new_Data, Blocks = prepro$Blocks, clust = clust, nstart = nstart,
    rho = rho, NameBlocks = prepro$NameBlocks, Itermax = Itermax,
    Graph_groups = Graph_groups, Graph_weights = Graph_weights, scale = FALSE,
    print_attempt = print_attempt
  )
  return(a)
}
