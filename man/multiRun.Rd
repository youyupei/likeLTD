\name{multiRun} 
\alias{multiRun} 
\title{likeLTD::multiRun} 
\usage{multiRun(hypothesis, nrun, \dots)}
\arguments{
  \item{hypothesis}{Hypothesis for which to run the optimisation}
  \item{nrun}{How many times to run the optimisation.}
  \item{\dots}{Any named parameter to be passed to optimisation.params() e.g. \code{fixed}}
}
\value{A list containing the standard deviation of the likelihood over nrun runs as well as the number of runs that completed. Additionally a set of parameters from the run that gave tha maximal likelihood; the starting parameters, the optimised parameters, the likelihood, and convergence information.} 
\description{Runs optimisation mulitple times}
\details{
  optimisation.params() selects random starting values for dropout and rcont. Different starting values may give slightly different maximal likelihoods (there may be multiple local maxima). multiRun() runs the optimisation with multiple starting values, and gives a measaure of how variable the results are.
}
\seealso{optimisation.params}
