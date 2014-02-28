\name{optimisation.params} 
\alias{optimisation.params} 
\title{likeLTD::optimisation.params} 
\usage{optimisation.params(hypothesis, verbose=TRUE, fixed=NULL,
                           logObjective=TRUE, logDegradation=TRUE,
                           arguments=NULL, zero=0, throwError=FALSE,
                           withPenalties=TRUE, objective=NULL, iterMax=NULL, 
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
  \item{objective}{Objective function produced from create.likelihood.vectors}
  \item{iterMax}{Number of iterations to run the optimisation for}
  \item{Whether to return likelihoods for every genotype combination, or 
  a likelihood summed over all genotypes after optimisation. Set to TRUE 
  for individual genotype likelihoods. This is used for 
  \code{get.likely.genotypes}.}
  \item{\dots}{Any named parameter to modify the hypothesis, e.g.
  \code{nUnknowns}}
}
\value{
  \item{par}{A vector of initial (nuisance) parameters}
  \item{fn}{The objective function}
  \item{lower}{Lower bounds for the parameters}
  \item{upper}{Upper bounds for the parameters}
  \item{control}{Control parameters for \code{stats::optim}}
  \item{method}{\dQuote{L-BFGS-B}}
  \item{hessian}{\code{FALSE}}
}
\description{Creates a list of parameters to use with \code{stats::optim}.}
\details{
  Starting from the hypothesis, it creates an list of arguments which can be
  applied to \code{stats::optim} to obtain the maximum (log-)likelihood of that
  hypothesis. By default, it provides the input for using the \dQuote{L-BFGS-B} of
  \code{optim}, including the upper and lower bounds of the nuisance parameter. 

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
