\name{initial.arguments} 
\alias{initial.arguments} 
\title{likeLTD::initial.arguments} 
\usage{initial.arguments(hypothesis, \dots)}
\arguments{
  \item{hypothesis}{Hypothesis from which objective function was obtained.}
  \item{\dots}{Any parameter to modify the hypothesis.}
}
\value{
  \item{rcont}{Vector of relative contributions. The relative contribution of
  the first known profile is set to 1 by default and thus not present in this
  vector.}
  \item{degradation}{Degradation for each contributor}
  \item{localAdjustment}{Adjustment for each locus}
  \item{tvedebrink}{Tvedebrink exponent}
  \item{dropin}{Dropin rate. Can be ignored if not dropin.}
  \item{dropout}{Dropout rates per replicate}
}
\description{Provides starting point values for optimization.}
\seealso{optimization.params}
