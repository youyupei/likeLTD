\name{optimisation.params} 
\alias{optimisation.params} 
\title{likeLTD::optimisation.params} 
\usage{optimisation.params(hypothesis, verbose=FALSE, fixed=NULL,
                           logObjective=TRUE, logDegradation=TRUE,
                           arguments=NULL, zero=0, throwError=FALSE,
                           withPenalties=TRUE, doLinkage=TRUE,objective=NULL, iterMax=75, 
                           likeMatrix=FALSE,\dots)}
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
                    prosecution likelihood by \eqn{IMP_L/IMP_U}, where \eqn{IMP_L} is the inverse match
                    probability with linkage taken into account, and \eqn{IMP_U} is the same
                    but without linkage taken into account.}
  \item{objective}{Objective function produced from create.likelihood.vectors}
  \item{iterMax}{Number of iterations to run the optimisation for}
  \item{likeMatrix}{Whether to return likelihoods for every genotype combination, or 
  a likelihood summed over all genotypes after optimisation. Set to TRUE 
  for individual genotype likelihoods. This is used for 
  \code{get.likely.genotypes}.}
  \item{\dots}{Any named parameter to modify the hypothesis, e.g.
  \code{nUnknowns}}
}
\value{
  \item{fn}{The objective function}
  \item{lower}{Lower bounds for the parameters}
  \item{upper}{Upper bounds for the parameters}
  \item{control}{Control parameters for \code{DEoptim::DEoptim}}
}
\description{Creates a list of parameters to use with \code{DEoptim::DEoptim}.}
\details{
  Starting from the hypothesis, it creates an list of arguments which can be
  applied to \code{DEoptim::DEoptim} to obtain the maximum (log-)likelihood of that
  hypothesis. 

  It accepts a number of customization:
    
  \itemize{
    \item The optimisation can be performed for the likelihood or the log of
    the likelihood. The latter is recommended.
    \item wether the \code{degradation} parameter should be inputs as \eqn{x}
    or  as an exponent \eqn{10^x}. The latter seems to be more numerically
    stable, likely because degradations (in first form) are factors of an
    exponent in any case.
    \item whether to keep some nuisance parameters fixed
  }

  In any case, the value returned can always be modified prior to calling
  \code{DEoptim::DEoptim}.
}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"hammer")

# File paths and case name for allele report
admin = pack.admin.input(
            cspFile = file.path(datapath, 'hammer-CSP.csv'),
            refFile = file.path(datapath, 'hammer-reference.csv'),
            caseName = "hammer",
	    kit= "SGMplus"
             )

# Enter arguments
args = list(
        nUnknowns = 1,
        doDropin = FALSE,
        ethnic = "EA1",
        adj = 1,
        fst = 0.02,
        relatedness = c(0,0)
        )

# Create hypotheses
hypP = do.call(prosecution.hypothesis, append(admin,args))
hypD = do.call(defence.hypothesis, append(admin,args))

# Get parameters for optimisation
paramsP = optimisation.params(hypP)
paramsD = optimisation.params(hypD)
}
}



