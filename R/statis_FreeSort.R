##' @title Performs the STATIS method on Free Sorting data
##'
##' @usage
##' statis_FreeSort(Data, NameSub=NULL, Graph_obj=TRUE, Graph_weights=TRUE)
##'
##' @description
##' STATIS method on Free Sorting data. A lot of supplementary informations are also computed
##'
##'
##' @param Data data frame or matrix. Corresponds to all variables that contain subjects results. Each column corresponds to a subject and gives the groups to which the products (rows) are assigned
##'
##' @param NameSub string vector. Name of each subject. Length must be equal to the number of clumn of the Data. If NULL, the names are S1,...Sm. Default: NULL
##'
##' @param Graph_obj logical. Show the graphical representation od the objects? Default: TRUE
##'
##' @param Graph_weights logical. Should the barplot of the weights be plotted? Default: TRUE
##'
##'
##'
##' @return a list with:
##' @return a list with:
##'         \itemize{
##'          \item RV: the RV matrix: a matrix with the RV coefficient between subjects
##'          \item compromise: a matrix which is the compromise of the subjects (akin to a weighted average)
##'          \item weights: the weights associated with the subjects to build the compromise
##'          \item lambda: the first eigenvalue of the RV matrix
##'          \item overall error : the error for the STATIS criterion
##'          \item error_by_conf: the error by configuration (STATIS criterion)
##'          \item rv_with_compromise: the RV coefficient of each subject with the compromise
##'          \item homogeneity: homogeneity of the subjects (in percentage)
##'          \item coord: the coordinates of each object
##'          \item eigenvalues: the eigenvalues of the svd decomposition
##'          \item inertia: the percentage of total variance explained by each axis
##'          \item error_by_obj: the error by object (STATIS criterion)
##'          \item scalefactors: the scaling factors of each subject
##'          \item proj_config: the projection of each object of each subject on the axes: presentation by subject
##'          \item proj_objects: the projection of each object of each subject on the axes: presentation by object
##'          }
##'
##'
##'
##'
##' @references
##' \itemize{
##' \item Lavit, C., Escoufier, Y., Sabatier, R., Traissac, P. (1994). The act (statis method). Computational 462 Statistics & Data Analysis, 18 (1), 97-119.\\
##' \item Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods.Application to sensometrics. Food Quality and Preference, in Press.
##' }
##'
##'
##' @keywords FreeSorting
##'
##' @examples
##'
##' data(choc)
##' res.sta=statis_FreeSort(choc)
##'
##' @seealso  \code{\link{preprocess_FreeSort}}, \code{\link{clustatis_FreeSort}}
##'
##' @export


statis_FreeSort <- function(Data, NameSub = NULL, Graph_obj = TRUE, Graph_weights = TRUE) {
  prepro <- preprocess_FreeSort(Data, NameSub = NameSub)

  a <- statis(
    Data = prepro$new_Data, Blocks = prepro$Blocks, NameBlocks = prepro$NameBlocks,
    scale = FALSE, Graph_obj = Graph_obj, Graph_weights = Graph_weights
  )
  return(a)
}
