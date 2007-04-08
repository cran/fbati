\name{fbati}
\alias{fbati}
\title{fbati}
\description{
  Family based test for gene-environment interaction for bi-allelic snps, command/line or GUI (provide no options to start gui, i.e. just type \code{fbati()} and press return).
}
\usage{
fbati( ped=NULL, phe=NULL,
       data=mergePhePed(ped,phe),
       marker=NULL, ## pairs...
       env,
       method="fbati",
       model="additive",
       iter=10000,
       seed=7,
       maxSib=3,
       fixNames=TRUE,
       debug=FALSE )
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped}. See \code{write.ped} in the pbatR R package for more details on the file format.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe}. See \code{write.phe} in the pbatR R package for more details on the file format.}
  \item{data}{a data.frame object containing required data, or formed from merging a pedigree and phenotype object together. The first columns of it must be as in a `ped' object, while the next can be in any order representing marker or phenotype information.}
  \item{marker}{Default is NULL for all markers. Otherwise, it can be the names of the marker (if you load in with read.ped, this should be without the '.a'/'.b' added to differentiate the two markers). If you are using more specialized loading routines, this represents the numbers of the columns where the markers are at. For example, 7:10 would mean that columns 7 and 8 represent one locus, and columns 9 and 10 represent another locus.}
  \item{env}{Character string representing the name of the environmental variable to use (a column header name of the 'data' parameter).}
  \item{method}{Currently only `fbati' is supported.}
  \item{model}{one of \code{"additive"}, \code{"dominant"}, or \code{"recessive".}}
  \item{iter}{The number of Monte-Carlo iterations to perform.}
  \item{seed}{The random seed, so consistent answers are maintained. See NOTE 1 for more details. NA/NULL disables this, but is not recommended.}
  \item{maxSib}{When nonzero, employs the following rules to minimize the number of strata, to improve the number of informative transmissions. When there are parents, a random affected child is chosen. When parents are missing, a random affected child with environmental exposure is chosen, and random genotyped siblings are chosen to maxSib total offspring (so 2 indicates a sibpair, 3 a sibtrio, etc.), and parents are treated as missing (even if there is one). See the 'strataReduce' routine for more details and examples.}
  \item{fixNames}{Just leave this to TRUE if creating from ped/phe objects (pops off the '.a' and '.b' added on to the names of the two markers that are added on when read in via the (f)read.ped(...) routine).}
  \item{debug}{Developer use only (extended output).}
}
\details{
  Returns a data.frame object with the results.  The columns entitled GX...X indicate the number of informative families in each strata for the given marker. If these columns do not show up, it indicates there was only one type of strata.

  The parents need not be in the dataset if they have completely missing genotypes (they will be inserted), but the snps must currently be bi-allelic (or you will get error messages).

  \code{fread.ped} and \code{fread.phe} are suggested to enforce loading the whole dataset.

  NOTE 1: The fbati test was developed for families with at least one affected, so if there is more than one affected individual per family, only a random affected one will be used, and a random unaffected to reduce strata, unless \code{strataFix} is disabled. This is done on a per marker basis, thus the seed is set before every marker to obtain reproducible results.

  NOTE 2: The data is converted into nuclear families. This is done by a
  call to `nuclifyMerged' to the passed in dataset to enforce this
  consistency.
}
\examples{
## Data is simulated according to the formula in the
##  paper (you can see it from the code).

## Set the seed so you get the same results
set.seed(13)

## Constants (you can vary these)
NUM_FAMILIES <- 500
AFREQ <- 0.1  ## Allele frequency
BG <- -0.25   ## main effect of gene
BE <- 0       ## main effect of environment
BGE <- 0.75   ## main gene-environment effect
ENV <- 0.2    ## environmental exposure frequency

## (but don't modify this one)
MAX_PROB <- exp( BG*2 + BE*1 + BGE*2*1 )

#####################################
## Create a random dataset (trios) ##
#####################################

## -- genotypes are set to missing for now,
##    everyone will be affected
ped <- as.ped( data.frame( pid = kronecker(1:NUM_FAMILIES,c(1,1,1)),
                           id = kronecker( rep(1,NUM_FAMILIES), c(1,2,3) ),
                           idfath = kronecker( rep(1,NUM_FAMILIES), c(0,0,1) ),
                           idmoth = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                           sex = rep(0,NUM_FAMILIES*3),
                           AffectionStatus = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                           m0.a = rep(0,NUM_FAMILIES*3),      ## missing for now
                           m0.b = rep(0,NUM_FAMILIES*3) ) )   ## missing for now
## -- envioronment not generated yet
phe <- as.phe( data.frame( pid = ped$pid,
                           id = ped$id,
                           env = rep(NA,NUM_FAMILIES*3) ) )   ## missing for now

## 50/50 chance of each parents alleles
mendelTransmission <- function( a, b ) {
  r <- rbinom( length(a), 1, 0.75 )
  return( a*r + b*(1-r) )
}

## Not the most efficient code, but it gets it done;
##  takes < 5 sec on pentium M 1.8Ghz
CUR_FAMILY <- 1
while( CUR_FAMILY<=NUM_FAMILIES ) {
  ## Indexing: start=father, (start+1)=mother, (start+2)=child
  start <- CUR_FAMILY*3-2

  ## Draw the parental genotypes from the population
  ped$m0.a[start:(start+1)] <- rbinom( 1, 1, AFREQ ) + 1
  ped$m0.b[start:(start+1)] <- rbinom( 1, 1, AFREQ ) + 1

  ## Draw the children's genotype from the parents
  ped$m0.a[start+2] <- mendelTransmission( ped$m0.a[start], ped$m0.b[start] )
  ped$m0.b[start+2] <- mendelTransmission( ped$m0.a[start+1], ped$m0.b[start+1] )

  ## Generate the environment
  phe$env[start+2] <- rbinom( 1, 1, ENV )

  ## and the affection status
  Xg <- as.integer(ped$m0.a[start+2]==2) + as.integer(ped$m0.b[start+2]==2)
  if( rbinom( 1, 1, exp( BG*Xg + BE*phe$env[start+2] + BGE*Xg*phe$env[start+2] ) / MAX_PROB ) == 1 )
    CUR_FAMILY <- CUR_FAMILY + 1
  ## otherwise it wasn't a valid drawn individual
}


##############
## Analysis ##
##############

## Print the first 4 families
print( head( ped, n=12 ) )
print( head( phe, n=12 ) )

## NOTE: We could have just put all of this info into a single dataframe otherwise,
##  that would look like just the results of this
data <- mergePhePed( ped, phe )
print( head( data, n=12 ) )

## And run the analysis on all the markers
fbati( ped=ped, phe=phe, env="env" )

## Or do it via the merged data.frame object
##  7 and 8 correspond to the marker columns
fbati( data=data, env="env", marker=c(7,8) )

## You may also want to up the number of Monte-Carlo
##  iterations from the default

## And we could also run a joint test instead
##  (see fbatj)
fbatj( ped=ped, phe=phe, env="env" )
fbatj( data=data, env="env", marker=c(7,8) )

\dontrun{
## This won't be run, but we could do this with the gui.
## It requires the data to be written to disk, so we do so:
write.ped( ped, "simulated" )
write.phe( phe, "simulated" )
## Then start the GUI -- specify the options as before,
##  but for the first two, navigate to the 'simulated.ped' and 'simulated.phe' files.
fbati()
}
}
\seealso{
  \code{\link{fbatj}},
}
\references{
Hoffmann, Thomas J.,  Lange, Christoph,  Vansteelandt, Stijn,  and Laird, Nan M.  Gene-Environment Interaction Test for Dichotomous Traits in Trios and Sibships. Submitted.

S. L. Lake and N. M. Laird. Tests of gene-environment interaction for case-parent triads with general environmental exposures. Ann Hum Genet, 68(Pt 1):55-64, Jan 2004.
}
\keyword{interface}
