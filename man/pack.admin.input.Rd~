\name{pack.admin.input} 
\alias{pack.admin.input} 
\title{likeLTD::pack.admin.input} 
\usage{pack.admin.input(cspFile, refFile, caseName='dummy',
                        databaseFile=NULL, kit=NULL, linkageFile=NULL,
                         outputPath=getwd())}
\arguments{
  \item{cspFile}{Mixed crime scene profile.}
  \item{refFile}{Reference profiles.}
  \item{caseName}{Name of the current case.}
  \item{databaseFile}{Path to the allele database.}
  \item{kit}{Database to use if \code{databaseFile} not specified.}
  \item{linkageFile}{Path to file with recombination rates for linked loci.}
  \item{outputPath}{Path where the output should be stored.}
}
\description{Packs and verifies administrative information. Only used in allele
report.}
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
}

