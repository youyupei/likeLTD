\name{Objective Functions} 
\alias{obfunc}
\title{likeLTD objective functions}
\usage{function(rcont, degradation, localAdjustment, tvedebrink, dropout,
                dropin, ...)}
\arguments{
  \item{rcont}{Relative contribution. Can be one of two lengths. See details.}
  \item{degradation}{Degradation parameters to determine likeliness of dropout}
  \item{localAdjustment}{A scalar (single locus objective) or a list of values
  giving the local adjustment for each locus.}
  \item{tvedebrink}{Tvedebrink exponent. Must be a scalar.}
  \item{dropout}{Dropout reate for each replicate}
  \item{dropin}{A scalar giving the dropin rate}
  \item{\dots}{Any other parameter. Mostly ignored by objective functions, but
  it makes it easier to pass on parameters to the penalty function}
}
\description{Objective functions giving the likelihoods computed by
\pkg{likeLTD}}
\details{
  The general workflow of \pkg{likeLTD} will generally result in the creation
  and maximization of one or more objective function. These functions will take
  the arguments listed above. They can also be accept any number of other
  arguments. This makes it easy both to extend \pkg{likeLTD} and to use
  objective functions with penalty functions.

  Some of the arguments are verified for size and/or value

  \itemize{
    \item{rcont}{Must be of length \eqn{k + u} or \eqn{k + u - 1}, where
    \eqn{k} is the number of contributors with known profiles, and \eqn{u} is
    the number of unknown contributors. Since the input is for \emph{relative}
    contributions, there are two possible sizes. In the first case, the
    contribution of the first known profile is kept to 1. In the second case,
    the contributions are used for all known and unknown contributors as given.
    It is for the user to know what she is doing. \code{rcont} should always be
    positive.}
    \item{degradation}{Should be of length \eqn{k + u}}
    \item{localAdjustment}{positive}
    \item{tvedebrink}{negative}
    \item{dropout}{\eqn{\in [0, 1]}{Between 0 and 1 included}}
    \item{dropin}{positive. Ignored if the scenario does not model dropin.}
  }
}
