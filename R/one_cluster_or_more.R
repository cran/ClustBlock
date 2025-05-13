.one_cluster_or_more <- function(Data, Blocks, nperm = 100, scale = FALSE, alpha = 0.05) {
  # initialisation
  n <- nrow(Data)
  nblo <- length(Blocks)
  J <- rep(1:nblo, times = Blocks)

  X <- scale(Data, center = TRUE, scale = scale) # X contains the centered (and scaled if necessary) data tables
  Wj <- array(0, dim = c(n, n, nblo)) # association matrices



  # globally standardization of each data matrix
  # Computation of association matrices
  for (j in 1:nblo) {
    Xj <- as.matrix(X[, J == j])
    Wj[, , j] <- tcrossprod(Xj)
    Wj[, , j] <- Wj[, , j] / sqrt(sum(diag(crossprod(Wj[, , j])))) # standardisation so that ||Wj||=1
  }

  vp <- .vp_rv_Wi(Wj)[2]


  lambda2_perm <- NULL
  for (h in 1:nperm)
  {
    Wj_per <- Wj
    for (k in 1:nblo)
    {
      valeur_pot <- Wj[, , k][upper.tri(Wj[, , k], diag = TRUE)]
      Wj_per[, , k][upper.tri(Wj_per[, , k], diag = TRUE)] <- sample(valeur_pot)
      for (j in 1:(n - 1))
      {
        for (i in (j + 1):n)
        {
          Wj_per[i, j, k] <- Wj_per[j, i, k]
        }
      }
    }
    lambda2_perm[h] <- .vp_rv_Wi(Wj_per)[2]
  }

  pval <- 0
  for (j in 1:nperm)
  {
    if (vp <= lambda2_perm[j]) {
      pval <- pval + 1
    }
  }
  pval <- pval / nperm
  decision <- pval > alpha
  return(list(decision = decision, pvalue = pval))
}
