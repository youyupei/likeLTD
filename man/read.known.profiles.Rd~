\name{read.known.profiles}
\alias{read.known.profiles}
\title{likeLTD::read.known.profiles}
\usage{read.known.profiles(path)}
\arguments{
  \item{path}{Path to a CSV file with known profiles. }
}
\value{
  A matrix where columns are loci and rows are individual profiles. Each
  element is a vector of two characters, and each character is the name of an
  allele represented in the individual. There should be two characters for each
  entry. 

  Additionnaly, the first column (named "queried") indicates whether that
  individiual is to be queried (TRUE) or is known to be in the CSP (FALSE).
}
\description{ Reads the known profiles from file. }
\seealso{read.csp.profile, read.unc.profiles}
\details{The input is a file in the CSV format (comma-separated values). 
It should have the following form:

  \tabular{llllccc}{
    Individual \tab known/queried \tab  \tab Locus \tab Locus \tab \dots \cr
    name \tab queried \tab "13,15" \tab "16,16" 
      \tab \dots\cr
    name \tab known \tab \tab "18, 19" \tab "16, 17.2"
      \tab \dots\cr
  }

Columns are separted by commas (not present in the table above). The first
line, containing the name of the columns need be present.  Each row is the
profile of a different individual. Names or identifiers for the individuals
are contained in the first column. The second column should contain the word
"queried" if the individual should be queried. The other columns are for the
differnt loci. "Locus" should be replaced with the name of the locus in the
first column. Each item consists of two comma separated values within qutoation
marks. The values are the name of the allles in that individual's profile for
that locus.
}
\examples{
######################
# Peak heights model #
######################
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"laboratory")

# File paths and case name for allele report
admin = pack.admin.input.peaks(
            peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
            refFile = file.path(datapath, 'laboratory-reference.csv'),
            caseName = "Laboratory",
            detectionThresh = 20
             )

# get reference profiles
read.known.profiles(admin$refFile)

##################
# Discrete model #
##################
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"hammer")

# File paths and case name for allele report
admin = pack.admin.input(
            cspFile = file.path(datapath, 'hammer-CSP.csv'),
            refFile = file.path(datapath, 'hammer-reference.csv'),
            caseName = "hammer",
	    kit= "SGMplus"
             )

# get reference profiles
read.known.profiles(admin$refFile)
}
