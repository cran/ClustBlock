## =============================================================================
##' @title Perform the CATATIS method on Just About Right data.
##'
##' @usage
##' catatis_jar(Data, nprod, nsub, levelsJAR=3, beta=0.1, Graph=TRUE, Graph_weights=TRUE,
##' Test_weights=FALSE, nperm=100)
##'
##' @description
##' CATATIS method adapted to JAR data.
##'
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
##' @param Graph logical. Show the graphical representation? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##' @param Test_weights logical. Should the the weights be tested? Default: FALSE
##'
##' @param nperm integer. Number of permutation for the weight tests. Default: 100
##'
##' @return a list with:
##'         \itemize{
##'          \item S: the S matrix: a matrix with the similarity coefficient among the subjects
##'          \item compromise: a matrix which is the compromise of the subjects (akin to a weighted average)
##'          \item weights: the weights associated with the subjects to build the compromise
##'          \item weights_tests: the weights tests results
##'          \item lambda:  the first eigenvalue of the S matrix
##'          \item overall error: the error for the CATATIS criterion
##'          \item error_by_sub: the error by subject (CATATIS criterion)
##'          \item error_by_prod: the error by product (CATATIS criterion)
##'          \item s_with_compromise: the similarity coefficient of each subject with the compromise
##'          \item homogeneity: homogeneity of the subjects (in percentage)
##'          \item CA: the results of correspondance analysis performed on the compromise dataset
##'          \item eigenvalues: the eigenvalues associated to the correspondance analysis
##'          \item inertia: the percentage of total variance explained by each axis of the CA
##'          \item scalefactors: the scaling factors of each subject
##'          \item nb_1: Can be ignored
##'          \item param: parameters called
##'          }
##'
##'
##'
##' @keywords JAR
##'
##' @references
##' Llobell, F., Vigneau, E. & Qannari, E. M. ((September 14, 2022). Multivariate data analysis and clustering of subjects in a Just about right task. Eurosense, Turku, Finland.
##'
##'
##'
##' @examples
##' data(cheese)
##' res.cat=catatis_jar(Data=cheese, nprod=8, nsub=72, levelsJAR=5)
##' summary(res.cat)
##' #plot(res.cat)
##'
##' @seealso   \code{\link{catatis}}, \code{\link{plot.catatis}}, \code{\link{summary.catatis}}, \code{\link{cluscata_jar}}, \code{\link{preprocess_JAR}}, \code{\link{cluscata_kmeans_jar}}
##'
##' @export

## =============================================================================
catatis_jar <- function(Data, nprod, nsub, levelsJAR = 3, beta = 0.1, Graph = TRUE, Graph_weights = TRUE, Test_weights = FALSE, nperm = 100) {
  # preprocessing
  prepro <- preprocess_JAR(Data, nprod, nsub, levelsJAR, beta)
  # catatis
  cat <- catatis(prepro$Datafinal, nsub, NameBlocks = prepro$NameSub, Graph = Graph, Graph_weights = Graph_weights, Test_weights = Test_weights, nperm = nperm)
  return(cat)
}
