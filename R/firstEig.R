.firstEig=function(mat, eps=1e-8)
{
  delta=1000
  n=nrow(mat)
  u=runif(n,-1000,1000)
 while(delta>eps)
 {
   unew=mat%*%u
   unew=unew/sqrt(sum(unew**2))
   delta=sqrt(sum((unew-u)**2))
   u=unew
 }
  lambda=(mat%*%u/u)[1]
  return(list(lambda=lambda,u=u))
}

