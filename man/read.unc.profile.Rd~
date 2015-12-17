\name{read.unc.profile}
\alias{read.unc.profile}
\title{likeLTD::read.unc.profile}
\usage{read.unc.profile(path)}
\arguments{
  \item{path}{Path to the crime scene profile. }
}
\value{
  A matrix where columns are loci and rows are replicates. Each element is a
  vector of character, and each character is the name of an allele \emph{which
  seems to be present but at a level too low for certainty} in the CSP for that
  loci and replicate.
}
\description{Reads the uncertain alleles from the Crime Scene Profile from
file. }
\details{The input is a CSV file. It is (or can be) the same file as for
\code{read.csv.profile}. In this case, the lines labeled "Uncertain" are
selected. 
}
\seealso{read.csp.profile, read.known.profiles}
\examples{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"hammer")

# File paths and case name for allele report
admin = pack.admin.input(
            cspFile = file.path(datapath, 'hammer-CSP.csv'),
            refFile = file.path(datapath, 'hammer-reference.csv'),
            caseName = "hammer",
	    kit= "SGMplus"
             )

# get uncertain profile
read.unc.profile(admin$cspFile)
}
