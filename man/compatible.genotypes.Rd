\name{compatible.genotypes} 
\alias{compatible.genotypes} 
\title{likeLTD::compatible.genotypes} 
\usage{compatible.genotypes(cspPresence, profPresence, alleleNames, nUnknowns,
                            dropin=FALSE, missingReps=NULL)}
\arguments{
  \item{cspPresence}{Crime Scene Profile for given locus. This should be a
        matrix where rows are replicas and contributors, columns are alleles
        for the given locus, and the input is \code{TRUE} or \code{FALSE} and
        indicates the presence of that particular allele.
  }
  \item{profPresence}{Known profiles for given locus. Same type of format as
        cspLocus.
  }
  \item{alleleNames}{Name of the alleles, e.g. columns of the previous two}
  \item{nUnknowns}{number of unknown contributors}
  \item{dropin}{\code{TRUE} if modelling drop-ins}
  \item{missingReps}{List of booleans indicating whether a replicate is missing
        in \code{cspPresence}
  }
}
\description{Genetics make-up of \code{nUnknown} unprofiled contributors for
given locus.}
\value{A m by (2n) matrix where the colums (grouped by twos) correspond to
       contributors, and each row is their potential contribution to the CSP.
}
