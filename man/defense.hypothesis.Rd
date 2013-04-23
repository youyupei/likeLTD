\name{defense.hypothesis}
\alias{defense.hypothesis}
\title{likeLTD::defense.hypothesis}
\description{Helper function to create the input for the defense.}
\usage{defense.hypothesis(mixedFile, refFile, ethnic='EA1', nUnknowns=0,
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
defense.}
\details{It loads the CSP, known profiles, and allele database from file. It
removes the queried individual from the known profiles. It increments the
number of unknown contributors by one (to make up for the queried individual).
The reference individual is set to "X" (first unprofiled contributor) by
default. 
} 
\seealso{prosecution.hypothesis}
