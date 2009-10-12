\name{fbatge}
\alias{fbatge}
\alias{fbatgeAll}
\title{fbatge}
\description{
  Family based test for gene-environment interaction utilizing arbitrary family structures and multiple affected offspring. This method is recommended over the \code{fbati} routine in most scenarios.

  If no arguments are passed, then a friendly graphical interface is presented to the user.

  fbatge [GxE test], fbatj (see fbatj help) [G,GxE test], fbatme (see fbatme help) [G test] generally have more options than fbatgeAll. fbatgeAll runs all three tests, and gives results of all of them, and so uses only the options that are common to all three functions.
}
\usage{
fbatge( ped=NULL, phe=NULL,
        env=NULL, cov=NULL,
        trait="AffectionStatus", geno=NULL,
        strategy="hybrid", model="additive" )
fbatgeAll( ped=NULL, phe=NULL, env=NULL, trait="AffectionStatus" )
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped} (non-symbolic). See \code{write.ped} in the pbatR R package for more details on the file format.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe} (non-symbolic). See \code{write.phe} in the pbatR R package for more details on the file format.}

  \item{env}{Environmental Exposure. Should be a string, corresponding to the name in the 'phe' object.}
  \item{cov}{Any covariates to adjust for (does not apply to RR method). Should be a vector of strings, corresponding to names in the 'phe' object.}
  \item{trait}{Dichotomous trait name. Should be either "AffectionStatus", corresponding to the affection status in the pedigree object, or a string in the phenotype object.}
  \item{geno}{Names of the genetic markers, from the 'ped' object. If NULL (default), then all genetic markers are tested.}
  \item{strategy}{One of 'hybrid' (recommended, most efficient, requires rare disease), 'RR' (relative risk model, generally for a rare disease), or 'CLR' (conditional logistic regression).}
  \item{model}{Either 'additive' for the additive genetic model, or 'codominant' for the codominant genetic model (indicator variables for the genotypes).}
}
\details{
  Implements the test as described in Hoffmann et. al. (please see References).

  NOTE: The allele frequency is simply based on the allele frequency in all genotyped individuals, and is not the best choice.
}
\examples{
  example( fbati ) ## See fbati, creates a dataset for us in 'phe' and 'ped'
  print( fbatge( ped=ped, phe=phe, env="env" ) )
  ## The results are very close to the FBAT-I function, which
  ##  we would expect for trios.
}
\references{
  Hoffmann, Thomas J., and Laird, Nan M. Combining Multiple Disease Models for a More Powerful Gene-Environment Interaction Test in Nuclear Families.
}
\keyword{interface}
