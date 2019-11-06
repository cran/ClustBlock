##=============================================================================


##' @title Perform a cluster analysis of free sorting data
##'
##' @description
##' Hierarchical clustering of free sorting data followed by a partitioning algorithm (consolidation). Each cluster of blocks is associated with a compromise
##' computed by the STATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##'clustatis_FreeSort(Data, NameSub=NULL, Noise_cluster=FALSE,Itermax=30,
##'                            Graph_dend=TRUE, Graph_bar=TRUE, printlevel=FALSE,
##'                            gpmax=min(6, ncol(Data)-1), Testonlyoneclust=TRUE,
##'                            alpha=0.05, nperm=50)
##'
##'
##' @param Data data frame or matrix. Corresponds to all variables that contain subjects results. Each column corresponds to a subject and gives the groups to which the products (rows) are assigned
##'
##' @param NameSub string vector. Name of each subject. Length must be equal to the number of clumn of the Data. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param Noise_cluster logical. Should a noise cluster be computed? Default: FALSE
##'
##' @param Itermax numerical. Maximum of iteration for the partitioning algorithm. Default: 30
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging be plotted? Default: FALSE
##'
##' @param printlevel logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default: min(6, ncol(Data)-1)
##'
##' @param Testonlyoneclust logical. Test if there is more than one cluster? Default: TRUE
##'
##' @param alpha numerical between 0 and 1. What is the threshold to test if there is more than one cluster? Default: 0.05
##'
##' @param nperm numerical. How many permutations are required to test if there is more than one cluster? Default: 50
##'
##'
##'
##'
##'
##'
##' @return Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
##'         \itemize{
##'          \item group: the clustering partition of subjects after consolidation. If Noise_cluster=TRUE, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold for the noise cluster
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item rv_with_compromise: RV coefficient of each block with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item comp_RV: RV coefficient between the compromises associated with the various clusters
##'          \item compromise: the W compromise of each cluster
##'          \item coord: the coordinates of objects of each cluster
##'          \item inertia: percentage of total variance explained by each axis for each cluster
##'          \item rv_all_cluster: the RV coefficient between each subject and each cluster compromise
##'          \item criterion: the CLUSTATIS criterion error
##'          \item param: parameters called in the consolidation
##'          \item type: parameter passed to other functions
##'          }
##'          There is also at the end of the list:
##'          \itemize{
##'          \item dend: The CLUSTATIS dendrogram
##'          \item cutree_k: the partition obtained by cutting the dendrogram for K clusters (before consolidation).
##'          \item overall_homogeneity_ng: percentage of overall homogeneity by number of clusters before consolidation (and after if there is no noise cluster)
##'          \item diff_crit_ng: variation of criterion when a merging is done before consolidation (and after if there is no noise cluster)
##'          \item test_one_cluster: decision and pvalue to know if there is more than one cluster
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
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
##'data(choc)
##'res.clu=clustatis_FreeSort(choc)
##'plot(res.clu, Graph_dend=FALSE)
##'summary(res.clu)
##'
##' @seealso   \code{\link{clustatis}}, \code{\link{preprocess_FreeSort}}, \code{\link{summary.clustatis}}, , \code{\link{plot.clustatis}}
##'
##' @export


##=============================================================================


clustatis_FreeSort=function(Data,NameSub=NULL, Noise_cluster=FALSE,Itermax=30,
                        Graph_dend=TRUE, Graph_bar=TRUE, printlevel=FALSE,
                        gpmax=min(6, ncol(Data)-1), Testonlyoneclust=TRUE,
                        alpha=0.05, nperm=50)
{

  prepro=preprocess_FreeSort(Data, NameSub = NameSub)

  a=clustatis(Data=prepro$new_Data,Blocks= prepro$Blocks,NameBlocks= prepro$NameBlocks,
              Noise_cluster=Noise_cluster, scale=FALSE, Itermax=Itermax,
              Graph_dend=Graph_dend, Graph_bar=Graph_bar, printlevel=printlevel, gpmax=gpmax,
              Testonlyoneclust=Testonlyoneclust, alpha=alpha, nperm=nperm)
  return(a)
}
