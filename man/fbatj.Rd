\name{fbatj}
\alias{fbatj}
\title{fbatj}
\description{
  Family based joint test for both the main effect of the gene and gene-environment interaction, command/line or GUI (provide no options to start gui, i.e. just type \code{fbati()} and press return).
}
\usage{
fbatj( ped=NULL, phe=NULL,
       data=mergePhePed(ped,phe),
       marker = NULL,
       trait = "AffectionStatus",
       env = "env",
       model = "additive",
       mode = "univariate",
       fixNames = TRUE,
       verbose = FALSE )
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped}. See \code{write.ped} in the pbatR R package for more details on the file format.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe}. See \code{write.phe} in the pbatR R package for more details on the file format.}
  \item{data}{a data.frame object containing required data, or formed from merging a pedigree and phenotype object together. The first columns of it must be as in a `ped' object, while the next can be in any order representing marker or phenotype information.}
  \item{marker}{Default is NULL for all markers. Otherwise, it can be the names of the marker (if you load in with read.ped, this should be without the '.a'/'.b' added to differentiate the two markers). If you are using more specialized loading routines, this represents the numbers of the columns where the markers are at. For example, 7:10 would mean that columns 7 and 8 represent one locus, and columns 9 and 10 represent another locus.}
  \item{trait}{Character string representing the name of the trait variable to use (a column header name of the 'data' parameter).}
  \item{env}{Character string representing the name of the environmental variable to use (a column header name of the 'data' parameter).}
  \item{model}{one of \code{"additive"}, \code{"dominant"}, or \code{"recessive".}}
  \item{mode}{"univariate" is the only one supported.}
  \item{fixNames}{Just leave this to TRUE if creating from ped/phe objects (pops off the '.a' and '.b' added on to the names of the two markers that are added on when read in via the (f)read.ped(...) routine).}
  \item{verbose}{Developer use only (extended output).}
}
\details{
  Returns a data.frame object with the results. Uses the empirical variance.

  For an example, see the fbati function.
}
\seealso{
  \code{\link{fbatj}}
}
\references{
K Lunetta, S V Faraone, J Biederman, and N M Laird. Family-based tests of association and linkage that use unaffected sibs, covariates, and interactions. Am J Hum Genet, 2000, 66, 605-614.
}