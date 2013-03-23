\name{determine.dropout}
\alias{determine.dropout}
\title{determine.dropout}
\usage{determine.dropout(knownProfiles, cspProfile)}
\arguments{
  \item{knownProfiles}{A matrix containing the profiles to examin.}
  \item{cspProfile}{A matrix containing the crime scene profile.}
}
\value{ A list with one element per individual. Each is True if that individual
 is subject to dropout. 
}
\description{An individual is subject to dropout if the individual's allele at
 one or more locus is not in the crime-scene profile.
}
\seealso{read.csp.profile, read.known.profiles}


