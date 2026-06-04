## =============================================================================
##' @title Perform a cluster analysis of variables
##'
##' @description
##' Perform cluster analysis of variables in context of quantitative,
##'  qualitative or mixed datasets with MixCluStatis.
##'
##' @usage
##' mixclustatis(Data, quanti=NULL,quali=NULL,Noise_cluster=FALSE,
##' Itermax=20, printlevel=FALSE, Graph_dend=TRUE, Graph_bar=TRUE,
##' gpmax=min(6, length(quali)+length(quanti)-1), rhoparam = NULL)
##'
##'
##' @param Data data frame or matrix. Correspond to all the data (variables are columns)
##'
##' @param quanti  numerical vector. The number of the columns containing quantitative variables.
##'
##' @param quali  numerical vector. The number of the columns containing qualitative variables.
##'
##' @param Noise_cluster logical. Should a noise cluster be computed? Default: FALSE
##'
##' @param Itermax numerical. Maximum of iteration for the partitioning algorithm. Default: 30
##'
##' @param printlevel logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Default: TRUE
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default:  min(6, number of variables -2)
##'
##' @param rhoparam numerical or vector. What is the threshold for the noise cluster? Between 0 and 1, high value can imply lot of vatriables set aside. If NULL, automatic threshold is computed. Can be different for each group (in this case, provide a vector)
##'
##'
##' @return Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
##'         \itemize{
##'          \item group: the clustering partition of variables after consolidation. If Noise_cluster=TRUE, some variables could be in the noise cluster ("K+1")
##'          \item rho: the threshold(s) for the noise cluster (computed or input parameter)
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item rv_with_compromise: RV coefficient of each variable with its cluster compromise
##'          \item weights: weight associated with each variable in its cluster
##'          \item comp_RV: RV coefficient between the compromises associated with the various clusters
##'          \item compromise: the W compromise of each cluster
##'          \item coord: the coordinates of objects of each cluster
##'          \item inertia: percentage of total variance explained by each axis for each cluster
##'          \item rv_all_cluster: the RV coefficient between each variable and each cluster compromise
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
##'          \item param: parameters called
##'          \item type: parameter passed to other functions
##'          }
##'
##'
##'
##' @keywords variables
##'
##' @references
##' Paper submitted: Llobell, F., Abdi, H., Eslami, A. (2026). Clustering of categorical and mixed data variables around latent variables.
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods. Application to sensometrics. Food Quality and Preference, in Press.\cr
##' Llobell, F., Vigneau, E., Qannari, E. M. (2019). Clustering datasets by means of CLUSTATIS with identification of atypical datasets. Application to sensometrics. Food Quality and Preference, 75, 97-104.
##' Llobell, F., & Qannari, E. M. (2020). CLUSTATIS: Cluster analysis of blocks of variables. Electronic Journal of Applied Statistical Analysis, 13(2).
##'
##'
##' @examples
##' data("wine", package = "FactoMineR")
##' res=mixclustatis(wine, quanti = 3:29, quali = 1:2)
##' summary(res)
##' plot(res, Graph_groups = FALSE, Graph_weights = TRUE)
##'
##' @seealso   \code{\link{clustatis}}, \code{\link{plot.clustatis}}, \code{\link{summary.clustatis}}
##'
##' @export


## =============================================================================


mixclustatis=function(Data, quanti=NULL,quali=NULL, Noise_cluster=FALSE,Itermax=20,
                       printlevel=FALSE, Graph_dend=TRUE, Graph_bar=TRUE,
                       gpmax=min(6, length(quali)+length(quanti)-1), rhoparam = NULL)
{

  if(is.null(quali)==TRUE & is.null(quanti)==TRUE)
  {
    quanti=1:ncol(Data)
  }


  p=length(quanti)+length(quali)
  n=nrow(Data)

  if(is.null(colnames(Data))) colnames(Data)=paste0("Y",1:ncol(Data))

  #qualitative data
  Yi=list()
  ind=0
  for (i in quali)
  {
    ind=ind+1
    Data[,i]=factor(Data[,i])
    Yi[[ind]]=as.matrix(tab.disjonctif(Data[,i]))
    a=apply(Yi[[ind]],2,mean)
    Yi[[ind]]=Yi[[ind]]%*%diag(1/sqrt(a))
    colnames(Yi[[ind]])=paste0("v",i,"m",1:ncol(Yi[[ind]]))
  }


  #quantitative data
  Xi=list()
  for (i in quanti)
  {
    Xi[[which(quanti==i)]]=Data[,i]
  }


  ind=0
  for (i in quali)
  {
    ind=ind+1
    Xi[[length(quanti)+ind]]=Yi[[ind]]
  }

  dat=NULL
  for (i in 1:p)
  {
    dat=cbind(dat,Xi[[i]])
  }
  rownames(dat)=rownames(Data)
  colum=sort(c(quanti,quali))
  Blocks=rep(0,p)
  Blocks[which(colum%in%quanti)]=1
  ind=0
  for (i in quali)
  {
    ind=ind+1
    Blocks[which(colum==quali[ind])]=nlevels(as.factor(Data[,i]))
  }

  colnam=colnames(Data[,sort(c(quanti,quali))])
  Data=dat
  NameBlocks =colnam

  a=clustatis(Data=dat,Blocks,NameBlocks =colnam,Noise_cluster=Noise_cluster,
              Itermax=Itermax,printlevel=printlevel, Graph_dend=Graph_dend,
              Graph_bar=Graph_bar, gpmax=gpmax, rhoparam = rhoparam,
              Testonlyoneclust = FALSE, scale=FALSE)
  return(a)
}

