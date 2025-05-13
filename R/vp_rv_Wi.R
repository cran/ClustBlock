.vp_rv_Wi <- function(Wi) {
  # Data size (n, all the var)
  # Blocks vector given the nb of var of each block

  # initialisation
  n <- dim(Wi)[1]
  nblo <- dim(Wi)[3]

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


  # eigenvalues of RV matrix
  ressvd <- svd(RV)
  lambda <- ressvd$d
  return(lambda[-nblo])
}
