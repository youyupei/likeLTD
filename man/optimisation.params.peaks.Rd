\name{optimisation.params.peaks} 
\alias{optimisation.params.peaks} 
\title{likeLTD::optimisation.params.peaks} 
\usage{optimisation.params.peaks(hypothesis, verbose=TRUE, fixed=NULL,
                           logObjective=TRUE, logDegradation=TRUE,
                           arguments=NULL, zero=1e-6, throwError=FALSE,
                           withPenalties=TRUE, doLinkage=TRUE,objective=NULL, iterMax=25, 
                           likeMatrix=FALSE,diagnose=FALSE,DEoptimStrategy=3,
			   searchPopFactor=1,DEoptimF=0.8,DEoptimC=NULL,
				maxDropin=250,\dots)}
\arguments{
  \item{hypothesis}{Hypothesis from which to perform maximization}
  \item{verbose}{Wether to print likelihood each and every time the objective
  function is called}
  \item{fixed}{Names of the parameters to keep fixed}
  \item{logObjective}{If \code{TRUE} (default), the objective function is the
  log10-likelihood.}
  \item{logDegradation}{If \code{TRUE} (default), the degradation parameters
  are entered as \eqn{10^x}}
  \item{arguments}{Initial parameters from which to start the maximization. If
  \code{NULL}, calls \code{initial.arguments}.}
  \item{zero}{Epsilon to indicate lower and upper bounds as \eqn{\alpha \pm
  \epsilon}{alpha +/- epsilon} that exclude the bound itself}
  \item{throwError}{If TRUE, throws an error if the result is infinite}
  \item{withPenalties}{If TRUE, then penalties are evaluated and used}
  \item{doLinkage}{Logical indicating whether or not to apply a correction for linked loci.
                    This correction is only applied when Q and X are assumed to be siblings
                    i.e. \code{hypothesis$relatedness=c(0.5,0.25)}. This multiplies the 
                    prosecution likelihood by \eqn{IMP_L/IMP_U}, where \eqn{IMP_L} is the 
		    inverse match probability with linkage taken into account, and \eqn{IMP_U} 
		    is the same but without linkage taken into account.}
  \item{objective}{Objective function produced from create.likelihood.vectors.peaks}
  \item{iterMax}{Number of iterations to run the optimisation for before checking for 
convergence within each step.}
  \item{likeMatrix}{Whether to return likelihoods for every genotype combination, or 
  a likelihood summed over all genotypes after optimisation. Set to TRUE 
  for individual genotype likelihoods. This is used for 
  \code{get.likely.genotypes.peaks}.}
  \item{diagnose}{Logical. If \code{TRUE} a list of peak heights, estimated mean peak heights
  and standard deviation of peak heights is returned by \code{result$fn}.}
  \item{DEoptimStrategy}{Parameter to control the strategy used by DEoptim. See 
	\code{DEoptim::DEoptim.control}.}
  \item{searchPopFactor}{The population size used by DEoptim is \code{searchPopFactor} 
	multiplied by the number of parameters to optimise over.}
  \item{DEoptimF}{Parameter to control differential weighting of DEoptim. See 
	\code{DEoptim::DEoptim.control}.}
  \item{DEoptimC}{Parameter to control the speed of crossover adaptation for DEoptim. See 
	\code{DEoptim::DEoptim.control}.}
  \item{maxDropin}{Upper bound for the dropin parameter to be used by \code{evaluate.peaks}.
	If you wish to explain a minor cnotributor as dropin this value may need to be greater
	than the default of 250.}
  \item{\dots}{Any named parameter to modify the hypothesis, e.g.
  \code{nUnknowns}}
}
\value{
  \item{fn}{The objective function}
  \item{lower}{Lower bounds for the parameters}
  \item{upper}{Upper bounds for the parameters}
  \item{control}{Control parameters for \code{DEoptim::DEoptim}}
}
\description{Creates a list of parameters to use with \code{DEoptim::DEoptim} for peak height data.}
\details{
  Starting from the hypothesis, it creates a list of arguments which can be
  applied to \code{DEoptim::DEoptim} to obtain the maximum (log-)likelihood of that
  hypothesis. 

  It accepts a number of customizations:
    
  \itemize{
    \item The optimisation can be performed for the likelihood or the log of
    the likelihood. The latter is recommended.
    \item whether the \code{degradation} parameter should be inputs as \eqn{x}
    or  as an exponent \eqn{10^x}. The latter seems to be more numerically
    stable, likely because degradations (in first form) are factors of an
    exponent in any case.
    \item whether to keep some nuisance parameters fixed.
  }

  In any case, the value returned can always be modified prior to calling
  \code{evaluate.peaks}.
}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"laboratory")

# File paths and case name for allele report
admin = pack.admin.input.peaks(
            peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
            refFile = file.path(datapath, 'laboratory-reference.csv'),
            caseName = "Laboratory",
            detectionThresh = 20
             )

# Enter arguments
args = list(
        nUnknowns = 1
        )

# Create hypotheses
hypP = do.call(prosecution.hypothesis.peaks, append(admin,args))
hypD = do.call(defence.hypothesis.peaks, append(admin,args))

# Get parameters for optimisation
paramsP = optimisation.params.peaks(hypP)
paramsD = optimisation.params.peaks(hypD)

# parameters without linkage adjustment
paramsP = optimisation.params.peaks(hypP,doLinkage=FALSE)
}
}


