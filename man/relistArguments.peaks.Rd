\name{relistArguments.peaks} 
\alias{relistArguments.peaks} 
\title{likeLTD::relistArguments.peaks} 
\usage{relistArguments.peaks(parameters, hypothesis, fixed=NULL, 
logDegradation=TRUE, arguments=NULL)}
\arguments{
  \item{parameters}{Vector of parameters}
  \item{hypothesis}{Hypothesis from which objective function was obtained.}
  \item{fixed}{Names of the arguments which were fixed during optimisation.}
  \item{logDegradation}{Whether degradation is logarithmic form.}
  \item{arguments}{Initial guess, if any, when starting minimization.}
}
\value{Input parameters as a list.}
\description{In practice, this function undoes the flattening needed to 
run DEoptim. Hence, it should take the same parameters that optimisation.params.peaks 
does. Takes a linear  vector of parameters, as passed to and returned from an 
optimisation method, and makes it a list. Transforms degradation
parameters back into normal (non-logarithmic) form. Adds fixed arguments back 
into the list.}
\seealso{optimisation.params.peaks, initial.arguments.peaks}



