## =============================================================================


##' @title Testing the difference in perception between two predetermined groups of subjects in a CATA experiment
##'
##' @description
##' Test adapted to CATA data to determine whether two predetermined groups of subjects have a different perception or not. For example, men and women.
##'
##' @usage
##'  simil_groups_cata(Data, groups, one=1, two=2, nperm=50, Graph=TRUE,
##'   alpha= 0.05, printl=FALSE)
##'
##'
##' @param Data data frame or matrix. Correspond to all the blocks of variables merged horizontally
##'
##' @param groups  categorical vector. The groups of each subject . The length must be the number of subjects.
##'
##' @param one string. Name of the group 1 in groups vector.
##'
##' @param two string. Name of the group 2 in groups vector.
##'
##' @param nperm numerical. How many permutations are required? Default: 50
##'
##' @param Graph logical. Should the CATATIS graph of each group be plotted? Default: TRUE
##'
##' @param alpha numerical between 0 and 1. What is the threshold of the test? Default: 0.05
##'
##' @param printl logical. Print the number of remaining permutations during the algorithm? Default: FALSE
##'
##'
##' @return a list with:
##'     \itemize{
##'     \item decision: the decision of the test
##'     \item pval: pvalue of the test
##'     }
##'
##'
##'
##' @keywords CATA
##'
##' @references
##' Llobell, F., Giacalone, D., Jaeger, S.R. & Qannari, E. M. (2021). CATA data: Are there differences in perception? JSM conference.\cr
##' Llobell, F., Giacalone, D., Jaeger, S.R. & Qannari, E. M. (2021). CATA data: Are there differences in perception? AgroStat conference.
##'
##'
##' @examples
##' \donttest{
##'  data(straw)
##'  groups=sample(1:2, 114, replace=TRUE)
##'  simil_groups_cata(straw, groups, one=1, two=2)
##' }
##' @export


## =============================================================================


simil_groups_cata <- function(Data, groups, one = 1, two = 2, nperm = 50, Graph = TRUE, alpha = 0.05, printl = FALSE) {
  nblo <- length(groups)
  nvar <- ncol(Data) / nblo
  Blocks <- rep(nvar, nblo)
  J <- rep(1:nblo, times = Blocks)
  n <- nrow(Data)

  if (length(unique(groups)) < 2) {
    stop("You need to have at least two groups in your groups vector")
  }

  # parapet for nblo
  if (as.integer(nvar) != nvar) {
    stop("number of columns modulo length of groups vector != 0")
  }

  # parapet for numerical Data
  for (i in 1:ncol(Data))
  {
    if (is.numeric(Data[, i]) == FALSE) {
      stop(paste("The data must be numeric (column", i, ")"))
    }
  }

  # parapet for binary Data
  if ((sum(Data == 0) + sum(Data == 1)) != (dim(Data)[1] * dim(Data)[2])) {
    stop("only binary Data is accepted (0 or 1)")
  }

  # parapet for number of objects
  if (n < 3) {
    stop("At least 3 products are required")
  }

  # parapet for number of blocks
  if (nblo < 4) {
    stop("At least 4 subjects are required")
  }

  # parapet for number of attributes
  if (nvar < 3) {
    stop("At least 3 attributes are required")
  }


  # no NA
  if (sum(is.na(Data)) > 0) {
    print("NA detected:")
    tabna <- which(is.na(Data), arr.ind = TRUE)
    print(tabna)
    stop(paste("NA are not accepted"))
  }


  ##### real value####
  # group 1
  cl1 <- which(groups == one)
  nblocl1 <- length(cl1)
  rescat1 <- catatis(Data[, J %in% cl1], nblocl1, Graph = FALSE, Graph_weights = FALSE)
  if (Graph == TRUE) {
    plot(rescat1, Graph_weights = FALSE, Graph_eig = FALSE, tit = paste(one))
  }

  # group 2
  cl2 <- which(groups == two)
  nblocl2 <- length(cl2)
  rescat2 <- catatis(Data[, J %in% cl2], nblocl2, Graph = FALSE, Graph_weights = FALSE)
  if (Graph == TRUE) {
    plot(rescat2, Graph_weights = FALSE, Graph_eig = FALSE, tit = paste(two))
  }
  Observed_value <- .s_between_comp(rescat1$compromise, rescat2$compromise)


  Sim_value <- NULL
  for (i in 1:nperm)
  {
    all <- sort(c(cl1, cl2))
    cl1perm <- sample(all, nblocl1)
    rescat1perm <- catatis(Data[, J %in% cl1perm], nblocl1, Graph = FALSE, Graph_weights = FALSE)
    # group 2
    cl2perm <- NULL
    for (j in 1:length(all))
    {
      if (all[j] %in% cl1perm == FALSE) {
        cl2perm <- c(cl2perm, all[j])
      }
    }
    rescat2perm <- catatis(Data[, J %in% cl2perm], nblocl2, Graph = FALSE, Graph_weights = FALSE)
    Sim_value[i] <- .s_between_comp(rescat1perm$compromise, rescat2perm$compromise)
    if (printl == TRUE) {
      print(i)
    }
  }

  pval <- sum(Sim_value < Observed_value) / nperm

  if (pval < alpha) {
    decision <- "The groups have difference in perception"
  } else {
    decision <- "The groups have no difference in perception"
  }

  return(list(decision = decision, pval = pval))
}
