## =============================================================================


##' @title Perform a cluster analysis of blocks of quantitative variables
##'
##' @description
##' Hierarchical clustering of quantitative Blocks followed by a partitioning algorithm (consolidation). Each cluster of blocks is associated with a compromise
##' computed by the STATIS method. Moreover, a noise cluster can be set up.
##'
##' @usage
##'  clustatis(Data,Blocks,NameBlocks=NULL,Noise_cluster=FALSE,
##'   Unique_threshold=TRUE,scale=FALSE,
##'   Itermax=30, Graph_dend=TRUE, Graph_bar=TRUE,
##'   printlevel=FALSE, gpmax=min(6, length(Blocks)-2),  rhoparam=NULL,
##'   Testonlyoneclust=FALSE, alpha=0.05, nperm=50)
##'
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param Blocks  numerical vector. The number of variables of each block. The sum must be equal to the number of columns of Data
##'
##' @param NameBlocks string vector. Name of each block. Length must be equal to the length of Blocks vector. If NULL, the names are B1,...Bm. Default: NULL
##'
##' @param Noise_cluster logical. Should a noise cluster be computed? Default: FALSE
##'
##' @param Unique_threshold logical. Use same rho for every cluster? Default: TRUE
##'
##' @param scale logical. Should the data variables be scaled? Default: FALSE
##'
##' @param Itermax numerical. Maximum of iteration for the partitioning algorithm. Default: 30
##'
##' @param Graph_dend logical. Should the dendrogram be plotted? Default: TRUE
##'
##' @param Graph_bar logical. Should the barplot of the difference of the criterion and the barplot of the overall homogeneity at each merging step of the hierarchical algorithm be plotted? Default: TRUE
##'
##' @param printlevel logical. Print the number of remaining levels during the hierarchical clustering algorithm? Default: FALSE
##'
##' @param gpmax logical. What is maximum number of clusters to consider? Default:  min(6, number of blocks -2)
##'
##' @param rhoparam numerical or vector. What is the threshold for the noise cluster? Between 0 and 1, high value can imply lot of blocks set aside. If NULL, automatic threshold is computed. Can be different for each group (in this case, provide a vector)
##'
##' @param Testonlyoneclust logical. Test if there is more than one cluster? Default: FALSE
##'
##' @param alpha numerical between 0 and 1. What is the threshold to test if there is more than one cluster? Default: 0.05
##'
##' @param nperm numerical. How many permutations are required to test if there is more than one cluster? Default: 50
##'
##'
##'
##' @return Each partitionK contains a list for each number of clusters of the partition, K=1 to gpmax with:
##'         \itemize{
##'          \item group: the clustering partition of datasets after consolidation. If Noise_cluster=TRUE, some blocks could be in the noise cluster ("K+1")
##'          \item rho: the threshold(s) for the noise cluster (computed or input parameter)
##'          \item homogeneity: homogeneity index (%) of each cluster and the overall homogeneity index (%) of the partition
##'          \item rv_with_compromise: RV coefficient of each block with its cluster compromise
##'          \item weights: weight associated with each block in its cluster
##'          \item comp_RV: RV coefficient between the compromises associated with the various clusters
##'          \item compromise: the W compromise of each cluster
##'          \item coord: the coordinates of objects of each cluster
##'          \item inertia: percentage of total variance explained by each axis for each cluster
##'          \item rv_all_cluster: the RV coefficient between each block and each cluster compromise
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
##'
##' @keywords quantitative
##'
##' @references
##' Llobell, F., Cariou, V., Vigneau, E., Labenne, A., & Qannari, E. M. (2018). Analysis and clustering of multiblock datasets by means of the STATIS and CLUSTATIS methods. Application to sensometrics. Food Quality and Preference, in Press.\cr
##' Llobell, F., Vigneau, E., Qannari, E. M. (2019). Clustering datasets by means of CLUSTATIS with identification of atypical datasets. Application to sensometrics. Food Quality and Preference, 75, 97-104.
##' Llobell, F., & Qannari, E. M. (2020). CLUSTATIS: Cluster analysis of blocks of variables. Electronic Journal of Applied Statistical Analysis, 13(2).
##'
##' @importFrom  stats as.dendrogram cor cutree runif
##'
##'
##' @examples
##'
##'  data(smoo)
##'  NameBlocks=paste0("S",1:24)
##'  cl=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks)
##'  #plot(cl, ngroups=3, Graph_dend=FALSE)
##'  summary(cl)
##'  #with noise cluster
##'  cl2=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks,
##'  Noise_cluster=TRUE, Graph_dend=FALSE, Graph_bar=FALSE)
##'  #with noise cluster and defined rho threshold
##'  cl3=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks,
##'  Noise_cluster=TRUE, Graph_dend=FALSE, Graph_bar=FALSE, rhoparam=0.5)
##'  #different Noise cluster thresholds
##'  cl4=clustatis(Data=smoo,Blocks=rep(2,24),NameBlocks = NameBlocks,
##'  Noise_cluster=TRUE, Graph_dend=FALSE, Graph_bar=FALSE, Unique_threshold= FALSE,
##'  rhoparam=c(0.6, 0.5,0.4))
##'
##' @seealso   \code{\link{plot.clustatis}}, \code{\link{summary.clustatis}} , \code{\link{clustatis_kmeans}}, \code{\link{statis}}
##'
##' @export


## =============================================================================





clustatis <- function(Data, Blocks, NameBlocks = NULL, Noise_cluster = FALSE,
                      Unique_threshold = TRUE, scale = FALSE,
                      Itermax = 30, Graph_dend = TRUE, Graph_bar = TRUE,
                      printlevel = FALSE, gpmax = min(6, length(Blocks) - 2), rhoparam = NULL,
                      Testonlyoneclust = FALSE, alpha = 0.05, nperm = 50) {
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
  if (nblo < 4) {
    stop("At least 4 blocks are required")
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

  # parapet for gpmax
  if (gpmax > (nblo - 1)) {
    stop(paste("gpmax > number of blocks-1"))
  }

  # parapet for rho
  if (Noise_cluster == TRUE && is.null(rhoparam) == FALSE) {
    if (Unique_threshold == TRUE) {
      if (length(rhoparam) > 1) {
        stop(paste("If Unique_threshold is TRUE, rhoparam should be NULL or a number between 0 and 1"))
      }
      if (rhoparam < 0 | rhoparam > 1) {
        stop(paste("If Unique_threshold is TRUE, rhoparam should be NULL or a number between 0 and 1"))
      }
    } else {
      if (length(rhoparam) > 1) {
        gpmax <- length(rhoparam)
      } else {
        stop(paste("If Unique_threshold is FALSE, rhoparam should be NULL or a vector of numbers between 0 and 1"))
      }
      for (i in 1:gpmax)
      {
        if (rhoparam[i] < 0 | rhoparam[i] > 1) {
          stop(paste("If Unique_threshold is FALSE, rhoparam should be NULL or a vector of numbers between 0 and 1"))
        }
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


  # centering and scaling if necessary
  Data <- scale(Data, center = TRUE, scale = scale)

  J <- rep(1:nblo, times = Blocks) # indicates which block each variable belongs to

  # parapet for constant configurations and compute of Wi

  Wi <- array(0, dim = c(nrow(Data), nrow(Data), nblo))
  for (i in 1:nblo)
  {
    Xi <- as.matrix(Data[, J == i])
    Wi[, , i] <- tcrossprod(Xi)
    nor <- sqrt(sum(diag(tcrossprod(Wi[, , i]))))
    if (nor == 0) {
      stop(paste("Configuration", i, "is constant"))
    }
    Wi[, , i] <- Wi[, , i] / nor # standardization
  }

  # only one cluster?
  if (Testonlyoneclust == TRUE) {
    testonecluster <- .one_cluster_or_more(Data, Blocks, nperm = nperm, scale = scale, alpha = alpha)
    if (testonecluster$decision == TRUE) {
      testonecluster$decision <- "Only one cluster can be considered"
    } else {
      testonecluster$decision <- "Clustering is necessary"
    }
  } else {
    testonecluster <- list(decision = "Untested", pvalue = "Untested")
  }


  # RV matrix:

  RV <- matrix(0, nblo, nblo)
  diag(RV) <- rep(1, nblo)
  if (nblo > 1) {
    for (i in 1:(nblo - 1)) {
      for (j in (i + 1):nblo) {
        RV[i, j] <- sum(diag(crossprod(Wi[, , i], Wi[, , j])))
        RV[j, i] <- RV[i, j]
      }
    }
  }



  # criterion Q when each data table forms a group itself: 0
  crit <- rep(0, nblo)
  Q_current <- 0
  # critS=rep(1,nblo)
  S_current <- nblo
  cc <- 1:nblo
  # the names of the clusters are 1 to 'number of clusters'
  code <- 1:nblo
  # the names of the clusters are taken from the HCA
  cvar <- 1:nblo
  # when a new cluster is formed from clusters i and j, it is named min(i,j))
  mergelist <- matrix(0, nblo - 1, 2)
  results <- matrix(0, nblo - 1, 5)
  colnames(results) <- c("merg1", "merg2", "new.clust", "agg.crit.hac", "clust.crit.hac")
  idgr <- matrix(0, nblo, nblo)
  idgr[nblo, ] <- 1:nblo
  ordr <- c()
  fusions <- -(1:nblo)
  hmerge <- matrix(0, nblo - 1, 2)
  ncluster <- nblo
  # for homogeneity
  quality <- NULL

  for (level in 1:(nblo - 1)) {
    ################
    # loops i et j
    ################

    # We search the two clusters that will be merged.

    deltamin <- 10^100
    for (i in 1:(ncluster - 1)) {
      for (j in (i + 1):(ncluster)) {
        # merging the i'th and j'th cluster:
        newcluster <- c(which(cc == i), which(cc == j))
        Wj <- list()
        for (tab in 1:length(newcluster)) {
          Wj[[tab]] <- Wi[, , newcluster[tab]]
        }

        newstatis <- .crit_statisWj_fast(Wj, newcluster, RV)
        deltacurrent <- newstatis - sum(crit[c(i, j)])
        # if deltacurrent is smaller than deltamin, the current information is saved:
        if (deltacurrent < deltamin) {
          deltamin <- deltacurrent
          c1 <- cvar[cc == i]
          c2 <- cvar[cc == j]
          merge1 <- c(c1[1], c2[1])
          c1 <- code[cc == i]
          c2 <- code[cc == j]
          merge <- c(c1[1], c2[1])
          cl_1 <- i
          cl_2 <- j
          statismerge <- newstatis
        }
      } # end of loop j
    } # end of loop i
    ncluster <- ncluster - 1

    ###############################
    # renewal of the parameters
    ###############################

    Q_current <- Q_current + deltamin
    mergelist[level, ] <- merge1
    results[level, 1:5] <- c(merge, nblo + level, deltamin, Q_current)
    crit[cl_1] <- statismerge
    crit <- crit[-cl_2]
    cc[which(cc == cl_2)] <- cl_1
    cc[which(cc > cl_2)] <- cc[which(cc > cl_2)] - 1
    cvar[cc == cl_1] <- min(merge1)
    code[cc == cl_1] <- nblo + level
    idgr[ncluster, ] <- cc
    indic <- c(fusions[merge1[1]], fusions[merge1[2]])
    hmerge[level, ] <- indic
    fusions[which(idgr[ncluster, ] == cl_1)] <- level
    ordr <- .order_var(ordr, which(idgr[ncluster, ] == cl_1))
    # #overall homogeneity
    quality <- c(quality, (nblo - Q_current) / nblo)

    # show the level
    if (printlevel == TRUE) {
      print(nblo - 1 - level)
    }
  } # end of level

  # Dendogram:
  resultscah <- list(labels = NameBlocks, height = results[, 4], merge = hmerge, order = ordr)
  mytot <- resultscah
  class(mytot) <- "hclust"
  mydendC <- as.dendrogram(mytot)

  # show the dendrogram
  if (Graph_dend == TRUE) {
    dev.new()
    cex <- 0.6
    par(cex = cex)
    par(mar = c(7, 4, 4, 2) + 0.1)
    plot(mydendC, type = "rectangle", main = "CLUSTATIS Dendrogram", axes = TRUE, cex = cex, ylab = "Height")
    par(cex = 1)
    par(mar = c(5, 4, 4, 2) + 0.1)
  }


  # show the criterion and homogeneity evolutions
  if (Graph_bar == TRUE) {
    dev.new()
    barplot(results[, 4][(nblo - gpmax):(nblo - 1)],
      xlab = "Nb clusters", ylab = "delta", main = "Variation of criterion",
      axisnames = TRUE, names.arg = paste((gpmax + 1):2, "->", (gpmax):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
    )

    dev.new()
    barplot(quality[(nblo - gpmax):(nblo - 1)] * 100,
      xlab = "Nb clusters", ylab = "Overall homogeneity (%)", main = "Overall homogeneity (%)",
      axisnames = TRUE, names.arg = paste((gpmax):1), las = 2, cex.names = 0.6, cex.main = 1.2, col = "blue"
    )
  }

  # number of clusters advised
  criter <- sort(results[, 5], decreasing = TRUE)
  H <- NULL
  for (k in 1:min(gpmax, nblo - 2))
  {
    H[k] <- (criter[k] / criter[k + 1] - 1) * (nblo - k - 1)
  }
  nbgroup_hart <- which.max(H[-(min(gpmax, nblo - 2))] - H[-1]) + 1
  cat(paste("Recommended number of clusters =", nbgroup_hart), "\n")





  ###### find the threshold(s) rho for the noise cluster and consolidation
  res.consol <- list()
  cutree_k <- list()
  rho <- NULL
  rho_k <- list()
  for (K in 1:gpmax)
  {
    coupe <- cutree(mytot, K)
    cutree_k[[K]] <- coupe
    if (Noise_cluster == TRUE) {
      if (is.null(rhoparam) == TRUE) {
        oldgroup <- coupe
        # compute the compromises
        Wk <- array(0, dim = c(nrow(Data), nrow(Data), K))
        for (i in 1:K)
        {
          cluster <- which(oldgroup == i)
          Wj <- list()
          for (tab in 1:length(cluster)) {
            Wj[[tab]] <- Wi[, , cluster[tab]]
          }
          statisk <- .crit_statisWj(Wj, cluster, RV)
          Wk[, , i] <- statisk$W
        }

        # compute the criterions
        cr <- list()
        cr2 <- list()
        for (i in 1:nblo)
        {
          a <- NULL
          for (k in 1:K)
          {
            W_k <- as.matrix(Wk[, , k])
            normW <- sum(diag(crossprod(W_k)))
            W_i <- Wi[, , i]
            a <- c(a, sum(diag(crossprod(W_i, W_k)))^2 / (normW))
          }
          cr[[i]] <- sqrt(a)
          if (K > 1) {
            cr2[[i]] <- sort(sqrt(a), decreasing = TRUE)[1:2]
          } else {
            cr2[[i]] <- sqrt(a)
          }
        }
        rho[K] <- mean(unlist(cr2))

        # multiple threshold
        if (Unique_threshold == FALSE) {
          grou=coupe

          mat_utile <- matrix(unlist(cr2), nrow = nblo, ncol = 2, byrow = TRUE)

          rho_k[[K]] <- 0
          for (k in 1:K)
          {
            rho_k[[K]] <- c(rho_k[[K]], (mean(mat_utile[which(grou == k), 1]) + mean(mat_utile[which(grou == k), 2])) / 2)
          }
          rho_k[[K]] <- rho_k[[K]][-1]
        }
      } else {
        if (Unique_threshold == TRUE) {
          rho[K] <- rhoparam
        } else {
          rho_k[[K]] <- rhoparam[1:K]
        }
      }
    } else {
      rho[K] <- 0
    }

    if (Unique_threshold == FALSE) {
      rho2consider <- unlist(rho_k[[K]])
    } else {
      rho2consider <- rho[K]
    }

    # consolidation

    res.consol[[K]] <- clustatis_kmeans(Data, Blocks, coupe,
      rho = rho2consider, NameBlocks = NameBlocks, Itermax = Itermax,
      Graph_groups = FALSE, Graph_weights = FALSE, scale = scale
    )
  }

  names(cutree_k) <- names(res.consol) <- paste0("partition", 1:gpmax)


  # after consolidation
  diff_crit_bef <- results[, 4][(nblo - gpmax + 1):(nblo - 1)]

  if (Noise_cluster == FALSE) {
    # after consolidation
    overall_after <- NULL
    crit_after <- NULL
    for (i in 1:gpmax)
    {
      overall_after[i] <- res.consol[[i]]$homogeneity[i + 1, 1]
      crit_after[i] <- res.consol[[i]]$criterion
    }
    diff_crit_after <- crit_after[-gpmax] - crit_after[-1]
    diff_crit_after <- sort(diff_crit_after)
    overall_after <- sort(overall_after, decreasing = TRUE)
    overall_homogeneity_ng <- rbind(quality[(nblo - gpmax):(nblo - 1)] * 100, overall_after)
    diff_crit_ng <- rbind(diff_crit_bef, diff_crit_after)
    rownames(diff_crit_ng) <- rownames(overall_homogeneity_ng) <- c("before consolidation", "after consolidation")
    colnames(diff_crit_ng) <- paste((gpmax):2, "->", (gpmax - 1):1)
    colnames(overall_homogeneity_ng) <- paste(gpmax:1, "cluster(s)")
  } else {
    overall_homogeneity_ng <- quality[(nblo - gpmax):(nblo - 1)] * 100
    diff_crit_ng <- diff_crit_bef
    names(diff_crit_ng) <- paste((gpmax):2, "->", (gpmax - 1):1)
    names(overall_homogeneity_ng) <- paste(gpmax:1, "cluster(s)")
  }



  # results
  res <- c(res.consol, list(
    dend = mydendC, cutree_k = cutree_k,
    overall_homogeneity_ng = round(overall_homogeneity_ng, 1),
    diff_crit_ng = round(diff_crit_ng, 2), test_one_cluster = testonecluster, param = list(
      nblo = nblo, ng = nbgroup_hart,
      Noise_cluster = Noise_cluster, gpmax = gpmax, n = n
    ), type = "H+C"
  ))
  class(res) <- "clustatis"

  return(res)
}
