% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consistency_cata.R
\name{consistency_cata}
\alias{consistency_cata}
\title{Test the consistency of each attribute in a CATA experiment}
\usage{
consistency_cata(Data,nblo, nperm=100, alpha=0.05, printAttrTest=FALSE)
}
\arguments{
\item{Data}{data frame or matrix. Correspond to all the blocks of variables merged horizontally}

\item{nblo}{numerical. Number of blocks (subjects).}

\item{nperm}{numerical. How many permutations are required? Default: 100}

\item{alpha}{numerical between 0 and 1. What is the threshold? Default: 0.05}

\item{printAttrTest}{logical. Print the number of remaining attributes to be tested? Default: FALSE}
}
\value{
a list with:
        \itemize{
         \item consist: the consistent attributes
         \item no_consist: the inconsistent attributes
         \item pval: pvalue for each test
         }
}
\description{
Permutation test on the agreement between subjects for each attribute in a CATA experiment
}
\examples{
\donttest{
 data(straw)
#with only 40 subjects
consistency_cata(Data=straw[,1:(16*40)], nblo=40)
#with all subjects
consistency_cata(Data=straw, nblo=114, printAttrTest=TRUE)
}

}
\references{
Llobell, F., Giacalone, D., Labenne, A.,  Qannari, E.M. (2019).	Assessment of the agreement and cluster analysis of the respondents in a CATA experiment.	Food Quality and Preference, 77, 184-190.
}
\seealso{
\code{\link{consistency_cata_panel}}, \code{\link{change_cata_format}}, \code{\link{change_cata_format2}}
}
\keyword{CATA}
\keyword{RATA}
