\name{defence.hypothesis.peaks}
\alias{defence.hypothesis.peaks}
\title{likeLTD::defence.hypothesis.peaks}
\description{Helper function to create the input for the defence using peak height data.}
\usage{defence.hypothesis.peaks(peaksFile, refFile, ethnic='NDU1', nUnknowns=0,
                        adj=1e0, fst=0.03, databaseFile=NULL, 
                        linkageFile=NULL,
                        relatedness=c(0,0), detectionThresh=30, doDropin=FALSE, 
                        doDoubleStutter=FALSE,doOverStutter=FALSE,
                        combineRare=TRUE, rareThreshold=0.05,
                        kit=NULL, ...)
}
\arguments{
  \item{peaksFile}{Path to the crime scene profile.}
  \item{refFile}{Path to the known profiles.}
  \item{ethnic}{Ethnicity, e.g. subgroup within the allele frequency database.}
  \item{nUnknowns}{Number of unknown contributors for which to perform
                   calculation.}
  \item{adj}{Allele frequency adjustment parameter.}
  \item{fst}{Allele fraction adjustment \eqn{F_{ST}}.}
  \item{databaseFile}{Path to the allele database. If NULL, then defaults to
                      the NGMSelect database provided with \code{likeLTD}. Must include
                      longest uninterrupted sequence (LUS) values for alleles.}
  \item{linkageFile}{Path to recombination rate information. If NULL then defaults
                    to the linkage file provided with \code{likeLTD}. Only used
                    if \code{relatedness}=\code{c(0.5,0.25)} i.e. when Q and X are 
                    siblings.}
  \item{relatedness}{IBD1 and IBD2 values to take into account close relatedness.}
  \item{detectionThresh}{Detection threshold for peaks.}
  \item{doDropin}{Whether or not to model drop-in. Note dropin is not currently possible
                    with the peak heights model.}
  \item{doDoubleStutter}{Logical. Whether or not to model double stutter.}
  \item{doOverStutter}{Logical. Whether or not to model over stutter.}                  
  \item{combineRare}{Whether or not to combine rare unobserved alleles into a single
                    allele.}
  \item{rareThreshold}{If \code{combineRare=TRUE}, this is the probability threshold
                        below which an allele is classed as rare, and therefore 
                        combined with other rare unobserved alleles.}
  \item{kit}{Parameter specifying which allele database supplied with likeLTD to use if 
            \code{linkageFile} is not specified. Possibilities are "DNA17-lus".}
  \item{...}{Other parameters to determine how to perform calculations.}
}
\value{A list of named input parameters, or hypothesis, suitable for the
defence.}
\details{It loads the CSP, known profiles, and allele database from file. It
removes the queried individual from the known profiles. It increments the
number of unknown contributors by one (to make up for the queried individual).
} 
\seealso{prosecution.hypothesis.peaks}


