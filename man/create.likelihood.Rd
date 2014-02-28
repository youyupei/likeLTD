\name{create.likelihood} 
\alias{create.likelihood} 
\title{likeLTD::create.likelihood} 
\usage{create.likelihood(hypothesis, addAttr=FALSE, \dots)}
\arguments{
  \item{hypothesis}{Hypothesis for which to create the objective function}
  \item{addAttr}{Whether to add attributes to the objective functions. There
  are two attributes: (1) \dQuote{hypothesis} referencing the input hypothesis,
  (ii) \sQuote{functions} containing the individual objective functions per
  Locus. The latter contain further attributes. }
  \item{\dots}{Any named parameter to modify the hypothesis, e.g.
  \code{nUnknowns}}
}
\value{A function} 
\description{Creates an objective function}
\details{
  The objective function is created from the hypothesis. Itself, it takes as
  arguments the nuisance parameters and, optionally, the parameters for the
  penalty function. This particular flavor of the objective returns the product
  of the likelihood and penalties across all loci.
}
\seealso{create.likelihood.vectors, create.likelihood.log, penalties, Objective Functions}
