.crit_statisWj_fast <- function(Wj, index, RVtot) {
  # Wj list

  n <- nrow(Wj[[1]])
  nblo <- length(Wj)

  RV <- RVtot[index, index]

  if (length(index) > 1) {
    ressvd <- svd(RV)
    lambda <- ressvd$d[1]
  } else {
    lambda <- 1
  }


  Q <- nblo - lambda

  return(Q)
}
