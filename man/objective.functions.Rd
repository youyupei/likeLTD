\name{Objective Functions} 
\alias{objective.function}
\title{likeLTD objective functions}
\description{Objective functions giving the likelihoods computed by
\pkg{likeLTD}.
\bold{These functions do not exist \emph{per se}. They are created by
\code{create.likelihood} and similar function.}}
\details{
  The general workflow of \pkg{likeLTD} will generally result in the creation
  and maximization of one or more objective function. These functions will take
  the arguments listed above. They can also be accept any number of other
  arguments. This makes it easy both to extend \pkg{likeLTD} and to use
  objective functions with penalty functions.

  Some of the arguments are verified for size and/or value. 
  
  The general form of these objective functions is:

  \code{function(rcont, degradation, localAdjustment, tvedebrink, dropout, dropin, \dots)}

  \itemize{
    \item{rcont}{Relative contribution. Must be of length \eqn{k + u} or \eqn{k
                 + u - 1}, where \eqn{k} is the number of contributors with
                 known profiles, and \eqn{u} is the number of unknown
                 contributors. Since the input is for \emph{relative}
                 contributions, there are two possible sizes. In the first
                 case, the contribution of a reference individual is kept to 1.
                 In general, the reference individual will be the queried
                 individual (if subjecto dropout and in prosecution hypothesis)
                 or the first individual subject to dropout.
                 In the second case, the contributions are used for all
                 known and unknown contributors as given.  It is for the user
                 to know what she is doing.  \code{rcont} should always be
                 positive.}
    \item{degradation}{Degradation parameters to determine likeliness of
                       dropout. Should be of length \eqn{k + u}.}
    \item{localAdjustment}{A scalar (single locus objective) or a list of values
                           giving the local adjustment for each locus. Must be
                           positive.}
    \item{tvedebrink}{Tvedebrink exponent. Must be a scalar and negative.}
    \item{dropout}{Dropout rate for each replicate. Must be between 0 and 1.}
    \item{dropin}{A scalar giving the dropin rate. Ignored in models without
                  dropin. Otherwise it should be positive (unless given in
                  exponential form).}
    \item{\dots }{Any other parameter. Mostly ignored by objective functions,
                 but it makes it easier to pass on parameters to the penalty
                 function}
  }
}
