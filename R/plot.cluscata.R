## =============================================================================

##' @title Displays the CLUSCATA graphs
##'
##' @usage
##' \method{plot}{cluscata}(x, ngroups=NULL, Graph_groups=TRUE, Graph_dend=TRUE,
##' Graph_bar=FALSE, Graph_weights=FALSE, axes=c(1,2), cex=1,
##' col.obj="blue", col.attr="red", ...)
##'
##' @description
##' This function plots dendrogram, variation of the merging criterion, weights and CATATIS map of each cluster
##'
##'
##' @param x object of class 'cluscata'.
##'
##' @param ngroups number of groups to consider. Ignored for cluscata_kmeans results. Default: recommended  number of clusters
##'
##' @param Graph_groups logical. Should each cluster compromise graphical representation be plotted? Default: TRUE
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Also available after consolidation if Noise_cluster=FALSE. Default: FALSE
##'
##' @param Graph_weights logical. Should the barplot of the weights in each cluster be plotted? Default: FALSE
##'
##' @param axes  numerical vector (length 2). Axes to be plotted. Default: c(1,2)
##'
##' @param cex numerical. Numeric character expansion factor; multiplied by par("cex") yields the final character size. NULL and NA are equivalent to 1.0.
##'
##' @param col.obj numerical or string. Color for the objects points. Default: "blue"
##'
##' @param col.attr numerical or string. Color for the attributes points. Default: "red"
##'
##' @param ... further arguments passed to or from other methods
##
##'
##'
##' @return the CLUSCATA graphs
##'
##'
##'
##'
##'
##'
##' @keywords CATA
##'
##' @examples
##' \donttest{
##'  data(straw)
##'  res=cluscata(Data=straw[,1:(16*40)], nblo=40)
##'  plot(res, ngroups=3, Graph_dend=FALSE)
##'  plot(res, ngroups=3, Graph_dend=FALSE,Graph_bar=FALSE, Graph_weights=FALSE, axes=c(1,3))
##' }
##' @seealso   \code{\link{cluscata}} , \code{\link{cluscata_kmeans}}
##'
##' @export


## =============================================================================





plot.cluscata <- function(x, ngroups = NULL, Graph_groups = TRUE, Graph_dend = TRUE,
                          Graph_bar = FALSE, Graph_weights = FALSE, axes = c(1, 2), cex = 1,
                          col.obj = "blue", col.attr = "red", ...) {
  res.cluscata <- x
  if (inherits(res.cluscata, "cluscata") == FALSE) {
    stop("The class of the object must be 'cluscata'")
  }

  n <- res.cluscata$param$n

  if (is.null(ngroups) | res.cluscata$type == "K") {
    ngroups <- res.cluscata$param$ng
  }

  if (res.cluscata$type == "H+C") {
    if (ngroups > res.cluscata$param$gpmax) {
      stop("ngroups>gpmax")
    }
  }

  # show the dendrogram
  if (Graph_dend == TRUE & res.cluscata$type == "H+C") {
    dev.new()
    par(cex = 0.6)
    par(mar = c(7, 4, 4, 2) + 0.1)
    plot(res.cluscata$dend, type = "rectangle", main = "CLUSCATA Dendrogram", axes = TRUE, cex = cex, ylab = "Height")
    par(cex = 1)
    par(mar = c(5, 4, 4, 2) + 0.1)
  }



  # show the criterion and homogeneity evolutions
  if (Graph_bar == TRUE & res.cluscata$type == "H+C") {
    nblo <- res.cluscata$param$nblo
    Noise_cluster <- res.cluscata$param$Noise_cluster
    gpmax <- res.cluscata$param$gpmax
    if (Noise_cluster == FALSE) {
      dif <- res.cluscata$diff_crit_ng[1, ]
      quality <- res.cluscata$overall_homogeneity_ng[1, ]
    } else {
      dif <- res.cluscata$diff_crit_ng
      quality <- res.cluscata$overall_homogeneity_ng
    }

    dev.new()
    barplot(dif,
      xlab = "Nb clusters", ylab = "delta", main = "Variation of criterion before consolidation",
      axisnames = TRUE, names.arg = paste((gpmax):2, "->", (gpmax - 1):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
    )

    dev.new()
    barplot(quality,
      xlab = "Nb clusters", ylab = "Overall homogeneity (%)", main = "Overall homogeneity before consolidation (%)",
      axisnames = TRUE, names.arg = paste((gpmax):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
    )

    if (Noise_cluster == FALSE) {
      dif <- res.cluscata$diff_crit_ng[2, ]
      quality <- res.cluscata$overall_homogeneity_ng[2, ]

      dev.new()
      barplot(dif,
        xlab = "Nb clusters", ylab = "delta", main = "Variation of criterion after consolidation",
        axisnames = TRUE, names.arg = paste((gpmax):2, "->", (gpmax - 1):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
      )

      dev.new()
      barplot(quality,
        xlab = "Nb clusters", ylab = "Overall homogeneity (%)", main = "Overall homogeneity after consolidation (%)",
        axisnames = TRUE, names.arg = paste((gpmax):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
      )
    }
  }

  if (res.cluscata$type == "H+C") {
    res.cluscata <- res.cluscata[[ngroups]]
  }

  if (Graph_groups == TRUE) {
    for (i in 1:ngroups)
    {
      un <- axes[1]
      deux <- axes[2]
      dev.new()
      options(ggrepel.max.overlaps = Inf)
      print(plot.CA(res.cluscata$CA[[i]],
        axes = axes, title = paste("Cluster", i),
        cex = cex, col.row = col.obj, col.col = col.attr
      ))
    }
  }

  if (Graph_weights == TRUE) {
    for (i in 1:ngroups)
    {
      dev.new()
      barplot(res.cluscata$weights[[i]])
      title(paste("Weight of each subject in the cluster", i))
    }
  }
}
