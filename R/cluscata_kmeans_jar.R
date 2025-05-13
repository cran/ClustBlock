## =============================================================================
##' @title Perform a cluster analysis of subjects in a JAR experiment
##' @description
##' Partitioning of subject from a JAR experiment. Each cluster is associated with a compromise
##' computed by the CATATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##' cluscata_kmeans_jar(Data, nprod, nsub, levelsJAR=3, beta=0.1, clust, nstart=100, rho=0,
##' Itermax=30, Graph_groups=TRUE, print_attempt=FALSE, Warnings=FALSE)
##'
##' @param Data data frame where the first column is the Assessors, the second is the products and all other columns the JAR attributes with numbers (1 to 3 or 1 to 5, see levelsJAR)
##'
##' @param nprod integer. Number of products.
##'
##' @param nsub integer. Number of subjects.
##'
##' @param levelsJAR integer. 3 or 5 levels. If 5, the data will be transformed in 3 levels.
##'
##' @param beta numerical. Parameter for agreement between JAR and other answers. Between 0 and 0.5.
##'
##' @param clust numerical vector or integer. Initial partition or number of starting partitions if integer. If numerical vector, the numbers must be 1,2,3,...,number of clusters
##'
##' @param nstart numerical. Number of starting partitions. Default: 100
##'
##' @param rho numerical or vector between 0 and 1. Threshold for the noise cluster. Default:0. If you want a different threshold for each cluster, you can provide a vector.
##'
##' @param Itermax numerical. Maximum of iterations by partitioning algorithm. Default: 30
##'
##' @param Graph_groups logical. Should each cluster compromise graphical representation be plotted? Default: TRUE
##'
##' @param print_attempt logical. Print the number of remaining attempts in multi-start case? Default: FALSE
##'
##' @param Warnings logical. Display warnings about the fact that none of the subjects in some clusters checked an attribute or product? Default: FALSE
##'
##' @return a list with:
##'         \itemize{
##'          \item group: the clustering partition. If rho>0, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold(s) for the noise cluster
##'          \item homogeneity: percentage of homogeneity of the subjects in each cluster and the overall homogeneity
##'          \item s_with_compromise: Similarity coefficient of each subject with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item compromise: The compromise of each cluster
##'          \item CA: The correspondance analysis results on each cluster compromise (coordinates, contributions...)
##'          \item inertia: percentage of total variance explained by each axis of the CA for each cluster
##'          \item s_all_cluster: the similarity coefficient between each subject and each cluster compromise
##'          \item param: parameters called
##'          \item criterion: the CLUSCATA criterion error
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##'
##' @keywords JAR
##'
##'
##' @references
##' Llobell, F., Vigneau, E. & Qannari, E. M. ((September 14, 2022). Multivariate data analysis and clustering of subjects in a Just about right task. Eurosense, Turku, Finland.
##'
##'
##' @examples
##' \donttest{
##' data(cheese)
##' res=cluscata_kmeans_jar(Data=cheese, nprod=8, nsub=72, levelsJAR=5, clust=4)
##' #plot(res)
##' summary(res)
##' }
##'
##' @seealso   \code{\link{plot.cluscata}}, \code{\link{summary.cluscata}} , \code{\link{catatis_jar}}, \code{\link{preprocess_JAR}}, \code{\link{cluscata_jar}}
##'
##' @export


## =============================================================================

cluscata_kmeans_jar <- function(Data, nprod, nsub, levelsJAR = 3, beta = 0.1, clust, nstart = 100, rho = 0, Itermax = 30,
                                Graph_groups = TRUE, print_attempt = FALSE, Warnings = FALSE) {
  # preprocessing
  prepro <- preprocess_JAR(Data, nprod, nsub, levelsJAR, beta)
  # catatis
  clu <- cluscata_kmeans(prepro$Datafinal, nsub,
    clust = clust, nstart = nstart, rho = rho, NameBlocks = prepro$NameSub, Itermax = Itermax,
    Graph_groups = Graph_groups, print_attempt = print_attempt, Warnings = Warnings
  )
  return(clu)
}
