## =============================================================================


##' @title Perform a cluster analysis of subjects from a RATA experiment
##'
##' @description
##' Hierarchical clustering of subjects (blocks) from a RATA experiment. Each cluster of blocks is associated with a compromise
##' computed by the CATATIS method. The hierarchical clustering is followed by a partitioning algorithm (consolidation).
##'
##' @usage
##' cluscata_rata(Data, nblo, NameBlocks=NULL, NameVar=NULL, Noise_cluster=FALSE,
##'         Unique_threshold =TRUE, Itermax=30, Graph_dend=TRUE,
##'         Graph_bar=TRUE, printlevel=FALSE,
##'         gpmax=min(6, nblo-2), rhoparam=NULL, Testonlyoneclust=FALSE, alpha=0.05,
##'         nperm=50, Warnings=FALSE)
##'
##' @param Data data frame or matrix where the blocks of binary variables are merged horizontally. If you have a different format, see \code{\link{change_cata_format}}
##'
##' @param nblo  numerical. Number of blocks (subjects).
##'
##' @param NameBlocks string vector. Name of each block (subject). Length must be equal to the number of blocks. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param NameVar string vector. Name of each variable (attribute, the same names for each subject). Length must be equal to the number of attributes. If NULL, the colnames of the first block are taken. Default: NULL
##'
##' @param Noise_cluster logical. Should a noise cluster be computed? Default: FALSE
##'
##' @param Unique_threshold logical. Use same rho for every cluster? Default: TRUE
##'
##' @param Itermax numerical. Maximum of iteration for the partitioning algorithm. Default:30
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Default: TRUE
##'
##' @param printlevel logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default: min(6, nblo-2)
##'
##' @param rhoparam numerical or vector. What is the threshold for the noise cluster? Between 0 and 1, high value can imply lot of blocks set aside. If NULL, automatic threshold is computed. Can be different for each group (in this case, provide a vector)
##'
##' @param Testonlyoneclust logical. Test if there is more than one cluster? Default: FALSE
##'
##' @param alpha numerical between 0 and 1. What is the threshold to test if there is more than one cluster? Default: 0.05
##'
##' @param nperm numerical. How many permutations are required to test if there is more than one cluster? Default: 50
##'
##' @param Warnings logical. Display warnings about the fact that none of the subjects in some clusters checked an attribute or product? Default: FALSE
##'
##'
##' @return Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
##'         \itemize{
##'          \item group: the clustering partition after consolidation. If Noise_cluster=TRUE, some subjects could be in the noise cluster ("K+1")
##'          \item rho: the threshold(s) for the noise cluster
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item s_with_compromise: similarity coefficient of each subject with its cluster compromise
##'          \item weights: weight associated with each subject in its cluster
##'          \item compromise: the compromise of each cluster
##'          \item CA: list. the correspondance analysis results on each cluster compromise (coordinates, contributions...)
##'          \item inertia: percentage of total variance explained by each axis of the CA for each cluster
##'          \item s_all_cluster: the similarity coefficient between each subject and each cluster compromise
##'          \item criterion: the CLUSCATA criterion error
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'          There is also at the end of the list:
##'          \itemize{
##'          \item dend: The CLUSCATA dendrogram
##'          \item cutree_k: the partition obtained by cutting the dendrogram in K clusters (before consolidation).
##'          \item overall_homogeneity_ng: percentage of overall homogeneity by number of clusters before consolidation (and after if there is no noise cluster)
##'          \item diff_crit_ng: variation of criterion when a merging is done before consolidation (and after if there is no noise cluster)
##'          \item test_one_cluster: decision and pvalue to know if there is more than one cluster
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##' @keywords RATA
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2019). A new approach for the analysis of data and the clustering of subjects in a CATA experiment. Food Quality and Preference, 72, 31-39.\cr
##' Llobell, F., Giacalone, D., Labenne, A.,  Qannari, E.M. (2019).	Assessment of the agreement and cluster analysis of the respondents in a CATA experiment.	Food Quality and Preference, 77, 184-190.
##' Llobell, F., Jaeger, S.R. (September 11, 2024). Consumer segmentation based on sensory product characterisations elicited by RATA questions? Eurosense conference, Dublin, Ireland.
##'
##' @examples
##'
##' #RATA data without session
##' data(RATAchoc)
##' Data=RATAchoc[1:108,2:16]
##' chang2=change_cata_format2(Data, nprod= 12, nattr= 13, nsub = 9, nsess = 1)
##' res.clus=cluscata_rata(Data= chang2$Datafinal, nblo = 9, NameBlocks =  chang2$NameSub)
##' summary(res.clus)
##' plot(res.clus)
##'
##'
##' @seealso   \code{\link{plot.cluscata}}, \code{\link{summary.cluscata}} , \code{\link{catatis_rata}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
##'
##' @export


## =============================================================================




cluscata_rata <- function(Data, nblo, NameBlocks = NULL, NameVar = NULL,
                          Noise_cluster = FALSE, Unique_threshold = TRUE, Itermax = 30,
                          Graph_dend = TRUE, Graph_bar = TRUE, printlevel = FALSE,
                          gpmax = min(6, nblo - 2), rhoparam = NULL,
                          Testonlyoneclust = FALSE, alpha = 0.05, nperm = 50, Warnings = FALSE) {
  # results
  res <- cluscata(
    Data, nblo, NameBlocks, NameVar, Noise_cluster, Unique_threshold,
    Itermax, Graph_dend, Graph_bar, printlevel, gpmax, rhoparam,
    Testonlyoneclust, alpha, nperm, Warnings
  )

  return(res)
}
