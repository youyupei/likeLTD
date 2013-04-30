\name{relistArguments} 
\alias{relistArguments} 
\title{likeLTD::relistArguments} 
\usage{relistArguments(arguments, hypothesis, logDegradation=TRUE)}
\arguments{
  \item{arguments}{Vector of arguments.}
  \item{hypothesis}{Hypothesis from which objective function was obtained.}
  \item{logDegradation}{Whether degradation is logarithmic form.}
}
\value{Input arguments as a list.}
\description{Takes a linear  vector of arguments, as passed to and returned
from an optimization method, and makes it a list. Transforms degradation
arguments back into normal (non-logarithmic) form. }
\seealso{initial.arguments, optimization.params}
