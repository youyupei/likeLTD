\name{load.allele.database} 
\alias{load.allele.database} 
\title{likeLTD::load.allele.database} 
\usage{load.allele.database(path=NULL, kit=NULL)}
\arguments{
  \item{path}{Path to the crime scene profile, or \code{NULL}.
              If \code{NULL}, then returns the allele database packaged with
              \pkg{likeLTD}.}
  \item{kit}{Only used if \code{path=NULL}. Specifies one of the allele database
            supplied with \code{likeLTD} to use. Possibilities are "DNA17",
            "Identifiler" and "SGMplus". Defaults to "DNA17" if both \code{path=NULL}
            and \code{kit=NULL}.}
}
\value{A table, as read from the file.}
\description{Reads allele database from file, or loads database provided with
\pkg{likeLTD}.}
\details{If not NULL, then the input should consist of a file in the following
format is using discrete model:

  \tabular{ccccccr}{
    Marker \tab Allele \tab  BP \tab EA1 \tab EA2 \tab EA3 \tab \dots \cr
    TH01 \tab 5   \tab 166 \tab   1  \tab  2  \tab   0 \cr
    TH01 \tab 6   \tab 170 \tab 212  \tab 44  \tab 118 \cr
    TH01 \tab 8   \tab 178 \tab  90  \tab 69  \tab  42 \cr
    TH01 \tab 8.3 \tab 181 \tab   0  \tab  0  \tab   1 \cr
    TH01 \tab 10  \tab 186 \tab   7  \tab  6  \tab   2 \cr
    TH01 \tab R   \tab  NA \tab   0  \tab  0  \tab   0 \cr
    vWA  \tab 13  \tab 163 \tab   0  \tab  5  \tab   2 \cr
    vWA  \tab 14  \tab 167 \tab  86  \tab 24  \tab  51 \cr
    \dots \cr
  }

The first line needs to be present. The first column indicated the locus, the
second the name of the allele, the third the fragment length. The next columns
are the frequencies for specific ethnic groups. There are no limits to the
number of ethnic groups included. If the name of an allele is
\code{R}, then it is ignored.

If using the peak heights model, there should be an extra column, "LUS", which
gives the mean longest uninterrupted sequence for that allele. 
}
\seealso{read.unc.profile, read.known.profiles}
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

# get allele database
load.allele.database(kit=admin$kit)
}
