\name{nuclify}
\alias{nuclify}
\alias{nuclifyMerged}
\alias{mergePhePed}
\title{Nuclify and Merge}
\description{
  \code{mergePhePed} merges a phenotype and pedigree object into a single data.frame object.

  \code{nuclifyMerged} chops a merged object into nuclear families of a dataset, generally a necessary preprocessing option for tests.

  \code{nuclify} chops instead a `ped' and `phe' object separately.
}
\usage{
mergePhePed( ped, phe )
nuclifyMerged( data )
nuclify( ped, phe )
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped}.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe}.}
  \item{data}{\code{data.frame} containing required data, or formed from merging a pedigree and phenotype object together. The first columns of it must be as in a `ped' object, while the next can be in any order representing marker or phenotype information.}
}
\details{
  \code{mergePhePed} and \code{nuclifyMerged} both return data.frame objects. \code{nuclify} returns a list that contains the `phe' object and the `ped' object with those respective names (see pbatR documentation, both objects extend data.frame objects, and can be used for the most part as if data.frame objects). When the data is nuclified, the parents of the nuclified families parents are lost.

  NOTE: nuclifyMerged will modify the pedigree id (pid) to be [100*(previous pid) + (nuclear family index)]. This should make it easy to observe the results of this call to your dataset.
}
\examples{
## Create some pedigree structure
##
##   100 --- 101
##        |
##       201---202
##           |
##     -------------
##     |   |   |   |
##    301 302 303 304
ped <- as.ped( data.frame( pid    = rep(1,8),
                           id     = c(100,101, 201,202, 301,302,303,304),
                           idfath = c(0,0,     100,0,   201,201,201,201),
                           idmoth = c(0,0,     101,0,   202,202,202,202),
                           sex    = c(1,2,     1,2,     2,2,2,2),
                           AffectionStatus = rep(0,8),
                           m0.a = rep(2,8),
                           m0.b = rep(2,8) ) )

## Which should chop up into
##   100 --- 101    201---202
##        |             |
##       201      -------------
##                |   |   |   |
##               301 302 303 304
nuclifyMerged( ped )

## NOTE: We could have merged the ped with a phe object,
##        via the 'mergePhePed' routine before running.
}
\keyword{interface}
