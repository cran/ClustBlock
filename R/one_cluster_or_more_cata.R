.one_cluster_or_more_cata <- function(Data, nblo, nperm = 100, alpha = 0.05) {
  n <- nrow(Data)
  p <- ncol(Data)
  nvar <- p / nblo
  Blocks <- rep(nvar, nblo)
  J <- rep(1:nblo, times = Blocks)
  vp <- .vp_s(Data, nblo)[2]
  lambda_per <- NULL
  for (i in 1:nperm)
  {
    Data_per <- Data
    for (j in 1:ncol(Data_per))
    {
      Data_per[, j] <- sample(Data[, j])
    }
    lambda_per[i] <- .vp_s(Data_per, nblo)[2]
  }

  pval <- 0
  for (i in 1:nperm)
  {
    if (vp <= lambda_per[i]) {
      pval <- pval + 1
    }
  }
  pval <- pval / nperm
  decision <- pval > alpha
  return(list(decision = decision, pvalue = pval))
}
