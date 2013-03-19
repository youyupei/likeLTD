\name{optimization.params} 
\alias{optimization.params} 
\title{likeLTD::optimization.params} 
\usage{optimization.params(scenario, verbose=TRUE, fixed=NULL,
                           logObjective=TRUE, logDegradation=TRUE,
                           arguments=NULL, zero=1e-8, \dots)}
\arguments{
  \item{scenario}{Scenario from which to perform maximixation}
  \item{verbose}{Wether to print likelihood each and every time the objective
  function is called}
  \item{fixed}{Names of the parameters to keep fixed}
  \item{logObjective}{If \code{TRUE} (default), the objective function is the
  log-likelihood.}
  \item{logDegradation}{If \code{TRUE} (default), the degradation parameters
  are entered as \eqn{10^x}}
  \item{arguments}{Initial parameters from which to start the maximization. If
  \code{NULL}, calls \code{initial.arguments}.}
  \item{zero}{Epsilon to indicate lower and upper bounds as \eqn{\alpha \pm
  \epsilon}{alpha +/- epsilon} that exclude the bound itself}
  \item{\dots}{Any named parameter to modify the scenario, e.g.
  \code{nUnknowns}}
}
\value{
  \item{par}{A vector of initial (nuisance) parametrs}
  \item{fn}{The objective function}
  \item{lower}{Lower bounds for the parameters}
  \item{upper}{Upper bounds for the parameters}
  \item{control}{Control parameters for \code{stats::optim}}
  \item{method}{\dQuote{L-BFGS-B}}
  \item{hessian}{\code{FALSE}}
}
\description{Creates a list of parameters to use with \code{stats::optim}.}
\details{
  Starting from the scenario, it creates an list of arguments which can be
  applied to \code{stats::optim} to obtain the maxmimum (log-)likelihood of that
  scenario. By default, it provides the input for using the \dQuote{L-BFGS-B} of
  \code{optim}, including the upper and lower bounds of the nuisance parameter. 

  It accepts a number of customization:
    
  \itemize{
    \item The optimization can be performed for the likelihood or the log of
    the likelihood. The latter is recommended.
    \item wether the \code{degradation} parameter should be inpute as \eqn{x}
    or  as an exponent \eqn{10^x}. The latter seems to be more numerically
    stable, likely because degradations (in first form) are factors of an
    expoenent in any case.
    \item whether to keep some nuisance parameters fixed
  }

  In any case, the value returned can always be modified prior to calling
  \code{stats::optim}.
}
