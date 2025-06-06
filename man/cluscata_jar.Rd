% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluscata_jar.R
\name{cluscata_jar}
\alias{cluscata_jar}
\title{Perform a cluster analysis of subjects in a JAR experiment.}
\usage{
cluscata_jar(Data, nprod, nsub, levelsJAR=3, beta=0.1,  Noise_cluster=FALSE,
        Unique_threshold=TRUE, Itermax=30, Graph_dend=TRUE, Graph_bar=TRUE,
         printlevel=FALSE, gpmax=min(6, nsub-2), rhoparam=NULL,
        Testonlyoneclust=FALSE, alpha=0.05, nperm=50, Warnings=FALSE)
}
\arguments{
\item{Data}{data frame where the first column is the Assessors, the second is the products and all other columns the JAR attributes with numbers (1 to 3 or 1 to 5, see levelsJAR)}

\item{nprod}{integer. Number of products.}

\item{nsub}{integer. Number of subjects.}

\item{levelsJAR}{integer. 3 or 5 levels. If 5, the data will be transformed in 3 levels.}

\item{beta}{numerical. Parameter for agreement between JAR and other answers. Between 0 and 0.5.}

\item{Noise_cluster}{logical. Should a noise cluster be computed? Default: FALSE}

\item{Unique_threshold}{logical. Use same rho for every cluster? Default: TRUE}

\item{Itermax}{numerical. Maximum of iteration for the partitioning algorithm. Default:30}

\item{Graph_dend}{logical. Should the dendrogram be plotted? Default: TRUE}

\item{Graph_bar}{logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Default: TRUE}

\item{printlevel}{logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE}

\item{gpmax}{logical. What is maximum number of clusters to consider? Default: min(6, nblo-2)}

\item{rhoparam}{numerical or vector. What is the threshold for the noise cluster? Between 0 and 1, high value can imply lot of blocks set aside. If NULL, automatic threshold is computed. Can be different for each group (in this case, provide a vector)}

\item{Testonlyoneclust}{logical. Test if there is more than one cluster? Default: FALSE}

\item{alpha}{numerical between 0 and 1. What is the threshold to test if there is more than one cluster? Default: 0.05}

\item{nperm}{numerical. How many permutations are required to test if there is more than one cluster? Default: 50}

\item{Warnings}{logical. Display warnings about the fact that none of the subjects in some clusters checked an attribute or product? Default: FALSE}
}
\value{
Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
        \itemize{
         \item group: the clustering partition after consolidation. If Noise_cluster=TRUE, some subjects could be in the noise cluster ("K+1")
         \item rho: the threshold(s) for the noise cluster
         \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
         \item s_with_compromise: similarity coefficient of each subject with its cluster compromise
         \item weights: weight associated with each subject in its cluster
         \item compromise: the compromise of each cluster
         \item CA: list. the correspondance analysis results on each cluster compromise (coordinates, contributions...)
         \item inertia: percentage of total variance explained by each axis of the CA for each cluster
         \item s_all_cluster: the similarity coefficient between each subject and each cluster compromise
         \item criterion: the CLUSCATA criterion error
         \item param: parameters called
         \item type: parameter passed to other functions
         }
         There is also at the end of the list:
         \itemize{
         \item dend: The CLUSCATA dendrogram
         \item cutree_k: the partition obtained by cutting the dendrogram in K clusters (before consolidation).
         \item overall_homogeneity_ng: percentage of overall homogeneity by number of clusters before consolidation (and after if there is no noise cluster)
         \item diff_crit_ng: variation of criterion when a merging is done before consolidation (and after if there is no noise cluster)
         \item test_one_cluster: decision and pvalue to know if there is more than one cluster
         \item param: parameters called
         \item type: parameter passed to other functions
         }
}
\description{
Hierarchical clustering of subjects from a JAR experiment. Each cluster of subjects is associated with a compromise
computed by the CATATIS method. The hierarchical clustering is followed by a partitioning algorithm (consolidation).
}
\examples{
\donttest{
data(cheese)
res=cluscata_jar(Data=cheese, nprod=8, nsub=72, levelsJAR=5)
#plot(res, ngroups=4, Graph_dend=FALSE)
summary(res, ngroups=4)
}

}
\references{
Llobell, F., Vigneau, E. & Qannari, E. M. ((September 14, 2022). Multivariate data analysis and clustering of subjects in a Just about right task. Eurosense, Turku, Finland.
}
\seealso{
\code{\link{plot.cluscata}}, \code{\link{summary.cluscata}} , \code{\link{catatis_jar}}, \code{\link{preprocess_JAR}}, \code{\link{cluscata_kmeans_jar}}
}
\keyword{JAR}
