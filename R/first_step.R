.first_step <- function(Data, nprod, nattr, nsub) {
  dat <- as.data.frame(Data)
  for (i in 0:(nprod - 1))
  {
    dat[seq(i + 1, nsub * nprod, nprod), ] <- Data[nsub * i + (1:nsub), ]
  }
  return(data = dat)
}
