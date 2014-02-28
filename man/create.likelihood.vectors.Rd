\name{create.likelihood.vectors} 
\alias{create.likelihood.vectors} 
\title{likeLTD::create.likelihood.vectors} 
\usage{create.likelihood.vectors(hypothesis, addAttr=FALSE, likeMatrix=FALSE, \dots)}
\arguments{
  \item{hypothesis}{Hypothesis for which to create the objective function}
  \item{addAttr}{Whether to add attributes to the objective functions. There
  are two attributes: (1) \dQuote{hypothesis} referencing the input hypothesis,
  (ii) \sQuote{functions} containing the individual objective functions per
  Locus. The latter contain further attributes. }
  \item{Whether to return likelihoods for every genotype combination, or 
  a likelihood summed over all genotypes after optimisation. Set to TRUE 
  for individual genotype likelihoods. This is used for 
  \code{get.likely.genotypes}.}
  \item{\dots}{Any named parameter to modify the hypothesis, e.g.
  \code{nUnknowns}}
}
\value{A function} 
\description{Creates an objective function}
\details{
  The objective function is created from the hypothesis. Itself, it takes as
  arguments the nuisance parameters and, optionally, the parameters for the
  penalty function. This particular flavor of the objective returns a list
  containing two items: (i) the likelihoods per locus, (ii) the penalties per
  locus.
}
\seealso{create.likelihood, create.likelihood.log, get.likely.genotypes, penalties, Objective Functions}
