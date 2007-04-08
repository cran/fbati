\name{strataReduce}
\alias{strataReduce}
\title{Strata Reduction}
\description{
  Reduces the number of strata in nuclear pedigrees for testing. For nuclear families with both parents, a random affected child is drawn. For nuclear families with at least one parent missing, a random affected and another random sib is used (parents ignored).
}
\usage{
strataReduce( data, envCol, m0, m1=m0+1, maxSib=3 )
}
\arguments{
  \item{data}{data.frame of a merged pedigree/phenotype (see mergePhePed(...)).}
  \item{envCol}{Integer representing environment column.}
  \item{m0}{Integer representing column of first marker.}
  \item{m1}{Integer representing column of second marker.}
  \item{maxSib}{Maximum number of sibs to use per family to reduce the number of strata.}
}
\examples{
## Function creates a family with the specified markers and statuses
## NOTE: affection is false/true, whereas it is coded 1/2 in the ped file
createFam <- function( pa=c(0,0), pb=c(0,0), ca, cb, caffected=rep(TRUE,length(ca)), env=1:length(ca) ) {
  ## pid, id, idfath, idmoth, sex, affection, m0a, m0b
  numC <- length(ca)
  return( data.frame( pid=rep(1,2+numC), id=1:(2+numC), idfath=c(0,0,rep(1,numC)), idmoth=c(0,0,rep(2,numC)), sex=c(2,1,rep(0,numC)), affection=c(0,0,as.integer(caffected)+1), m0.a=c(pa,ca), m0.b=c(pb,cb), env=c(NA,NA,env) ) )
}
## Function tests/exemplifies the strataReduce(...) routine
srFam <- function( ... ) {
  data <- createFam( ... )
  data2 <- strataReduce( data=data, envCol=9, m0=7, maxSib=2 )
  cat( "Original data:\n" )
  print( data )
  cat( "Reduced stratification data:\n" )
  print( data2 )
}

## Basic sib test
srFam( ca=c(1,1,2), cb=c(1,2,2) )

## Basic trio test
srFam( ca=c(1,1,2), cb=c(1,2,2), pa=c(1,1), pb=c(2,2) )

## a fairly comprehensive test here
## The affected should always be one of the first three,
##  the unaffected could be one the first eight
for( i in 1:10 )
  srFam( ca=c(1:8,0,0), cb=c(1:8,0,0), pa=c(1,1), caffected=c(rep(TRUE,6),rep(FALSE,4)), env=c(1:3,rep(NA,7)) )

## Now just to make sure, a full pedigree, rather than just one family
data <- createFam( ca=1:2, cb=1:2 )
for( i in 2:10 )
  data <- rbind( data, createFam( ca=1:2, cb=1:2 ) )
  cat( "Original data (full pedigree):\n" )
  print( data )
  cat( "Reduced stratification data (full pedigree), maxSib=3\n" )
  print( strataReduce( data=data, envCol=9, m0=7 ) )
}
\keyword{interface}
