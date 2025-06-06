## =============================================================================

##' @title Change format of CATA datasets to perform CATATIS or CLUSCATA function
##'
##' @usage
##' change_cata_format(Data, nprod, nattr, nsub, format=1, NameProds=NULL, NameAttr=NULL)
##'
##' @description
##' CATATIS and CLUSCATA operate on data where the blocksvariables are merged horizontally.
##' If you have a different format, you can use this function to change the format.
##' Format=1 is for data merged vertically with the dataset of the first subject, then the second,... with products in same order
##' Format=2 is for data merged vertically with the dataset for the first product, then the second... with subjects in same order
##'
##' Unlike change_cata_format2, you don't need to specify products and subjects, just make sure they are in the right order.
##'
##' @param Data data frame or matrix. Correspond to your data
##'
##' @param nprod integer. Number of products
##'
##' @param nattr integer. Number of attributes
##'
##' @param nsub integer. Number of subjects.
##'
##' @param format integer (1 or 2). See the description
##'
##' @param NameProds string vector with the names of the products (length must be nprod)
##'
##' @param NameAttr string vector with the names of attributes (length must be nattr)
##'
##'
##'
##' @return The arranged data for CATATIS and CLUSCATA function
##'
##'
##'
##' @keywords CATA
##'
##' @seealso   \code{\link{catatis}}, \code{\link{cluscata}}, \code{\link{change_cata_format2}}
##'
##' @export






change_cata_format <- function(Data, nprod, nattr, nsub, format = 1, NameProds = NULL, NameAttr = NULL) {
  if (format == 1) {
    return(.second_step(Data, nprod, nattr, nsub, NameProds, NameAttr))
  } else if (format == 2) {
    data <- .first_step(Data, nprod, nattr, nsub)
    return(.second_step(data, nprod, nattr, nsub, NameProds, NameAttr))
  }
}
