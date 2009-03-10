\name{fbatc}
\alias{fbatc}
\title{fbatc}
\description{
  Family based test for a group of markers conditional on another group of markers (typically conditional on a single marker). To start the graphical interface, provide no options, i.e. type \code{fbatc()} and press return.
}
\usage{
fbatc( ped=NULL, phe=NULL, data=mergePhePed( ped, phe),
       trait="AffectionStatus", traitType="auto",
       markerAnalyze=NULL, markerCondition=NULL,
       offset=NULL,
       tempPrefix="temp",
       MAXITER=1000, TOL=sqrt(.Machine$double.eps),
       verbose=FALSE )
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped} (non-symbolic). See \code{write.ped} in the pbatR R package for more details on the file format.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe} (non-symbolic). See \code{write.phe} in the pbatR R package for more details on the file format.}
  \item{data}{a data.frame object containing required data, or formed from merging a pedigree and phenotype object together. The first columns of it must be as in a `ped' object, while the next can be in any order representing marker or phenotype information.}

  \item{trait}{Trait to be analyzed. Defaults to AffectionStatus.}
  \item{traitType}{``auto'',``binary'', or ``continuous'': if set to ``auto'', then ``binary'' will be chosen if there is only two levels of outcome, otherwise ``continuous''.}

  \item{markerAnalyze}{Names of markers to analyze (without .a, e.g.).}
  \item{markerCondition}{Names of markers to condition on. If none are specified, then each marker will be conditioned on in turn.}
  \item{offset}{If set to NULL (i.e. left unset, the default) then for a continuous trait this is estimated by the trait mean.}

  \item{tempPrefix}{Temporary prefix to use for output files. These are safe to delete later.}

  \item{MAXITER}{Maximum iterations before giving up on convergence for the nuisance parameters.}
  \item{TOL}{Relative tolerance for convergence for the nuisance parameters.}

  \item{verbose}{For debug.}
}
\details{
  Implements the test as described in Hoffmann et. al. (please see References).

  The results returned are a data.frame object. The column `pvalue' and `rank' are the pvalue and rank of the empirical covariance matrix of the model-based test (dichotomous or normal). The column `pvalueR' and `rankR' are the pvalue and rank of the robust test. The model-based test has considerable more power over the robust test, but must assume a disease model. Please see Hoffmann et. al. for more details.

  This requires that FBAT be installed. If it is not, then the routine will attempt to automatically install it when given permission to do so by the user.
}
\examples{
\dontrun{
  set.seed(13)

  ## We simulate NO LD HERE, and a completely random trait!
  ## Data here is only to show how to use the function

  ###################
  ## IGNORE START: ##
  ###################

  ##  You can safely ignore how the data is generated,
  ##  and just see how to use it afterward.
  NUM_FAMILIES <- 500
  AFREQ <- c(0.2,0.2)
  ped <- as.ped( data.frame( pid = kronecker(1:NUM_FAMILIES,c(1,1,1)),
                            id = kronecker( rep(1,NUM_FAMILIES), c(1,2,3) ),
                            idfath = kronecker( rep(1,NUM_FAMILIES), c(0,0,1) ),
                            idmoth = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                            sex = rep(0,NUM_FAMILIES*3),
                            AffectionStatus = kronecker( rep(1,NUM_FAMILIES), c(0,0,2) ),
                            m0.a = rep(0,NUM_FAMILIES*3),      ## missing for now
                            m0.b = rep(0,NUM_FAMILIES*3),
                            m1.a = rep(0,NUM_FAMILIES*3),
                            m1.b = rep(0,NUM_FAMILIES*3) ) )
  CUR_FAMILY <- 1
  while( CUR_FAMILY<=NUM_FAMILIES ) {
    ## Indexing: start=father, (start+1)=mother, (start+2)=child
    start <- CUR_FAMILY*3-2

    ## Draw the parental genotypes from the population
    ped$m0.a[start:(start+1)] <- rbinom( 1, 1, AFREQ[1] ) + 1
    ped$m0.b[start:(start+1)] <- rbinom( 1, 1, AFREQ[1] ) + 1
    ped$m1.a[start:(start+1)] <- rbinom( 1, 1, AFREQ[2] ) + 1
    ped$m1.b[start:(start+1)] <- rbinom( 1, 1, AFREQ[2] ) + 1

    ## Draw the children's genotype from the parents
    ma <- rbinom( 1, 1, 0.5 )
    mb <- rbinom( 1, 1, 0.5 )
    if( rbinom( 1, 1, 0.5 ) == 0 ) {
      ped$m0.a[start+2] <- ped$m0.a[start]
      ped$m1.a[start+2] <- ped$m1.a[start]
    }else{
      ped$m0.a[start+2] <- ped$m0.b[start]
      ped$m1.a[start+2] <- ped$m1.b[start]
    }
    if( rbinom( 1, 1, 0.5 ) == 0 ) {
      ped$m0.b[start+2] <- ped$m0.a[start+1]
      ped$m1.b[start+2] <- ped$m1.a[start+1]
    }else{
      ped$m0.b[start+2] <- ped$m0.b[start+1]
      ped$m1.b[start+2] <- ped$m1.b[start+1]
    }

    CUR_FAMILY <- CUR_FAMILY + 1
  }

  ## Create a completely random phenotype as well
  phe <- as.phe( data.frame( pid=ped$pid, id=ped$id, qtl=rnorm(nrow(ped)) ) )

  ################
  ## IGNORE END ##
  ################

  ## Look at the first part of the pedigree
  print( head( ped ) )
  ## Look at the first part of the phenotype
  print( head( phe ) )

  ## Binary trait
  ## -- fbatc default trait is AffectionStatus
  ## -- fbatc default trait type is 'auto'
  ## - Test marker m1 conditional on m0
  print(  fbatc( ped=ped, markerAnalyze="m1", markerCondition="m0" )  )
  ## - Do the test the other way around, m0 conditional on m1
  print(  fbatc( ped=ped, markerAnalyze="m0", markerCondition="m1" )  )
  ## - Otherwise, we could have done this in one step;
  ##    if markerCondition is not specified, each member
  ##     of markerAnalyze is used.
  print(  fbatc( ped=ped, markerAnalyze=c("m0","m1") )  )

  ## QTL
  print(  fbatc( ped=ped, phe=phe, trait="qtl", markerAnalyze="m1", markerCondition="m0" )  )
  print(  fbatc( ped=ped, phe=phe, trait="qtl", markerAnalyze="m0", markerCondition="m1" )  )

  ## Additionally, we could write out the data that we
  ##  generated to disk so that we can then use it.
  write.ped( "simulated", ped ) ## to simulated.ped
  write.phe( "simulated", phe ) ## to simulated.phe
}
}
\references{
  Hoffmann, Thomas J. and Laird, Nan M. Parsing the Effects of Individual SNPs in Candidate Genes in Families. Submitted.
}
\keyword{interface}
