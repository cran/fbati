\name{fbatcStrategyStep}
\alias{fbatcStrategyStep}
\alias{fbatcStrategyStepUp}
\alias{fbatcStrategyStepDown}
\alias{fbatcStrategyStepLatex}
\alias{print.fbatcSStep}
\title{FBAT-C Stepwise Strategy}
\description{
	Apply the FBAT-C test in a stepwise fashion using \code{fbatcStrategyStep} (which does forward selection with \code{fbatcStrategyStepUp}, followed by backwards selection with \code{fbatcStrategyStepDown}) and get the results ready for publication with \code{fbatcStrategyStepLatex}.
}
\usage{
fbatcStrategyStepUp( ped, phe, markers=ped.markerNames(ped), trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE )
fbatcStrategyStepDown( ped, phe,  markers=ped.markerNames(ped), markersChosen=ped.markerNames(ped), markersChosenR=markersChosen, trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE )
fbatcStrategyStep( ped, phe, markers=ped.markerNames(ped), trait="trait", traitType="auto", alphaMMarker=0.05, alphaStep=alphaMMarker, sortByCorrelation=TRUE, tempPrefix="temp_", sim=FALSE, debug=FALSE )
fbatcStrategyStepLatex( res, digits=4, ffile="", preamble=FALSE, build=preamble, pdf="" )
\method{print}{fbatcSStep}(x,...)
}
\arguments{
  \item{ped}{Object from \code{(f)read.ped} or \code{as.ped}. See \code{write.ped} in the pbatR R package for more details on the file format.}
  \item{phe}{Object from \code{(f)read.phe} or \code{as.phe}. See \code{write.phe} in the pbatR R package for more details on the file format.}
  \item{markers}{Names of the markers to analyze.}
  \item{trait}{Name of the trait to analyze. Can be dichotomous or continuous.}
  \item{traitType}{"auto","dichotomous", or "continuous". If "auto" (the default), then "dichotomous" will be set if there are only two levels of the phenotype.}
	\item{alphaMMarker}{Alpha value for the multimarker test.}
	\item{alphaStep}{Alpha value used in the stepwise procedure.}
	\item{sortByCorrelation}{Whether to sort the markers by putting those in highest correlation closest to each other.}
	\item{tempPrefix}{The prefix to use for some intermittent files. Changing this is only necessary when you want to run this routine in parallel when each process shares the same disk.}
	\item{sim}{Developer use only.}
	  
	\item{res}{Result of 'fbatcStrategyStep' routine.}
	\item{digits}{Number of significant digits to display.}
	\item{ffile}{If set to a filename, then the output is redirected to that file instead of the standard output.}
	\item{preamble}{Whether to produce a latex file that can be compiled, or only the code for the chart.}
	\item{build}{Whether to run pdflatex on the file (requires \code{preamble=TRUE}), pdflatex must be in your path (generally true in linux, but not in Windows).}
	\item{pdf}{Name of the pdf viewer executable, if you also want to open the compiled file immediately. Note that in this case, you may not be able to return to the R session until you close this window.}
	
	\item{markersChosen}{In the step-down approach, the markers to start with for the model-based approach.}
	\item{markersChosenR}{In the step-down approach, the markers to start with for the model-free approach.}
	
	\item{debug}{Developer use only.}

  \item{x}{Result of fbatcStrategyStep, fbatcStrategyStepUp, fbatcStrategyStepDown.}
  \item{...}{Extra arguments.}
}
\details{
	fbatcStrategy returns a list with the following components.
	
	\code{mmarkerPvalue}: p-value of the multi-marker test on those markers (Rakovski et. al 2008).
	
	\code{correlation}: correlation matrix of the markers
	
	\code{univariate}: univariate results
	
	\code{step}: (model-based test) list of components \code{pvalue} (ith pvalue of the conditional test of markersAnalyze[i] on all markersCondition), \code{numInf} (number of informative families in the ith test), \code{markersAnalyze}, and \code{markersCondition}
	
	\code{markersChosen}: (model-based test) results from applying step-up strategy
	
	\code{stepR}, \code{markersChosenR}: (model-free test) results similar to \code{step} and \code{markersChosen}.
}
\keyword{interface}
