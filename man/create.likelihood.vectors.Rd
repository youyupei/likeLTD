\name{create.likelihood.vectors} 
\alias{create.likelihood.vectors} 
\title{likeLTD::create.likelihood.vectors} 
\usage{create.likelihood.vectors(scenario, addAttr=FALSE, \dots)}
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
  penalty function. This particular flavor of the objective returns a list
  containing two items: (i) the likelihoods per locus, (ii) the penalties per
  locus.
}
\seealso{create.likelyhood, create.likelihood.log, penalties, Objective Functions}
