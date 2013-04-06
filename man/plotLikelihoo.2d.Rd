\name{plotLikelihood.2d} 
\alias{plotLikelihood.2d} 
\title{likeLTD::plotLikelihood.2d} 
\usage{plotLikelihood.2d(hypothesis, which=c(1, 2), large=100, N=20,
                         arguments=NULL, x=NULL, y=NULL, logObjective=TRUE,
                         logDegradation=TRUE, contours=list())}
\arguments{
  \item{hypothesis}{Hypothesis from which to derive the likelihood}
  \item{which}{Indices of the arguments for which to plot likelihood}
  \item{large}{Number with which to replace Inf in upper and lower bounds. Can
               be a 2-tuple.}
  \item{N}{Number of points per axis. Can be a 2-tuple.}
  \item{arguments}{Arguments to the likelihood. Defines values for dimensions
                   that are not plotted. Defaults to
                   \code{initial.arguments(...)}.} 
  \item{x}{Range over which to plot abcissa} 
  \item{y}{Range over which to plot ordinate} 
  \item{logObjective}{Whether to plot log-likelihood}
  \item{logDegradation}{Whether to input degradation as \eqn{10^x}}
  \item{contours}{Arguments to stat_contour}
}
\description{Plots two dimensional graph of likelihood for two input
             directions.}
