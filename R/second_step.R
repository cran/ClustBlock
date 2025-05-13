.second_step <- function(Data, nprod, nattr, nsub, NameProds = NULL, NameAttr = NULL) {
  if (is.null(NameProds)) NameProds <- paste0("P", 1:nprod)
  if (is.null(NameAttr)) NameAttr <- paste0("A", 1:nattr)
  if (length(NameProds) != nprod) {
    stop("Length of NameProds must be equal to nprod")
  }
  if (length(NameAttr) != nattr) {
    stop("Length of NameAttr must be equal to nattr")
  }
  catadata <- as.data.frame(Data)
  dat <- matrix(0, nrow = nprod, ncol = nattr * nsub)
  dat <- as.data.frame(dat)
  for (i in 0:(nsub - 1))
  {
    dat[, (i * nattr + 1):(i * nattr + nattr)] <- catadata[(i * nprod + 1):(i * nprod + nprod), ]
  }
  colnames(dat) <- rep(NameAttr, nsub)
  rownames(dat) <- NameProds
  return(dat)
}
