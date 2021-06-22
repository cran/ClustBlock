.s_between_comp=function(comp1, comp2)
{
  num=sum(diag(tcrossprod(comp1, comp2)))
  den=sqrt(sum(diag(tcrossprod(comp1)))*sum(diag(tcrossprod(comp2))))
  return(num/den)
}
