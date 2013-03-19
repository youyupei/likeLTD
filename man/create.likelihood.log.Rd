\name{create.likelihood.log} 
\alias{create.likelihood.log} 
\title{likeLTD::create.likelihood.log} 
\usage{create.likelihood.log(scenario, addAttr=FALSE, \dots)}
\arguments{
  \item{scenario}{Scenario for which to create the objective function}
  \item{addAttr}{Whether to add attributes to the objective functions. There
  are two attributes: (1) \dQuote{scenario} referencing the input scenario,
  (ii) \sQuote{functions} containing the individual objective functions per
  Locus. The latter contain further attributes. }
  \item{\dots}{Any named parameter to modify the scenario, e.g.
  \code{nUnknowns}}
}
\value{A function} 
\description{Creates an objective function}
\details{
  The objective function is created from the scenario. Itself, it takes as
  arguments the nuisance parameters and, optionally, the parameters for the
  penalty function. This particular flavor of the objective returns the sum of
  the log-likelihood and the sum of log-penalties across loci.
}
\seealso{create.likelyhood.vectors, create.likelihood, penalties, Objective Functions}
