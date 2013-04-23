\name{prosecution.hypothesis}
\alias{prosecution.hypothesis}
\title{likeLTD::prosecution.hypothesis}
\usage{prosecution.hypothesis(mixedFile, refFile, ethnic='EA1', nUnknowns=0,
                              adj=1e0, fst=0.02, databaseFile=NULL, ...)
}
\arguments{
  \item{mixedFile}{Path to the crime scene profile.}
  \item{refFile}{Path to the known profiles.}
  \item{ethnic}{Ethnicity, e.g. subgroup within the allele frequency database.}
  \item{nUnknowns}{Number of unknown contributors for which to perform
                   calculation.}
  \item{adj}{Allele frequency adjustment parameter.}
  \item{fst}{Allele fraction adjustment \eqn{F_{ST}}}
  \item{databaseFile}{Path to the allele database. If NULL, then defaults to
                      the database provided with \code{likeLTD}}
  \item{...}{Other parameters to determine how to perform calculations, e.g.
             \code{doDropin}}
}
\value{A list of named input parameters, or hypothesis, suitable for the
prosecution.}
\description{Helper function to create the input for the prosecution.}
\details{It loads the CSP, known profiles, and allele database from file. It
makes sure the relatedness factor is zero. It sets the reference individual to
the first queried individual, unless it is not subject to dropout. In that
case, the reference individual is the first individual with dropout. } 
\seealso{defense.hypothesis}
