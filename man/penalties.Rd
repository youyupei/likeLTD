\name{penalties} 
\alias{penalties} 
\title{likeLTD::penalties} 
\usage{penalties(rcont, degradation, localAdjustment, tvedebrink, dropout,
                 dropin, localAdjPenalty=50, dropinPenalty=2,
                 degradationPenalty=50, bemn=-4.35, besd=0.38, ...)}
\arguments{
  \item{rcont}{Ignored.}
  \item{degradation}{Degradation parameters}
  \item{localAdjustment}{Local adjustment for each locus}
  \item{tvedebrink}{Tvedebrink exponent}
  \item{dropout}{Ignored}
  \item{dropin}{Dropin rate}
  \item{localAdjPenalty}{Penalty parameter for the local adjustments}
  \item{dropinPenalty}{Penalty parameter for the dropin rate}
  \item{degradationPenalty}{Penalty parameter for the degradation parameters}
  \item{bemn}{Mean of the normal distribution used to penalize degradation}
  \item{besd}{Standard deviation of the normal distribution used to penalize
  degradation}
  \item{\dots}{Ignored}
}
\description{Returns the penalty for the current arguments}
\value{An array of penalties per locus}
\details{
  The penalties are applied if and only if the relevant arguments
  (localAdjustment, dropin, degradation, tvedebrink) are provided. 
  The penalties are as follows:

  \itemize{
    \item{dropin:}{\eqn{e^{-d*p}} where \eqn{d} is the dropin rate
    and \eqn{p} the associated penalty. The values is normalized to one locus.}
    \item{degradation:}{\eqn{e^{-p\sum d}} where \eqn{d} are the degradation
    values and \eqn{p} is the associated penalty}
    \item{tvedebrink:}{\code{dnorm(t, bemn, besd)} where \code{t} is the
    Tvedebrink exponent, \code{dnorm} is the density of the normal distribution
    with mean \code{bemn} and standard deviation \code{besd}} 
    \item{localAdjustment:}{\code{dgamma(a, p, p)} where \code{a} is the local
    adjustment, \code{dgamma} is the density of the Gamma distribution with
    \code{p} its shape and rate}
  }

  Some of these penalties are meant to be applied simultaneously across all
  loci. Since we want penalties \emph{per} locus, a normalization
  \eqn{p^{\frac{1}{n}}}{p^(1/n)} is applied, where \eqn{p} is the penalty and
  \eqn{n} the number of loci.
}
\seealso{create.likelyhood.vectors, create.likelihood.log, create.likelyhood,
         Objective Functions}
