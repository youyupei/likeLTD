\name{allele.report} 
\alias{allele.report} 
\title{likeLTD::allele.report} 
\usage{allele.report(admin,file=NULL)}
\arguments{
  \item{admin}{List containing administration data, as packed by pack.admin.input()}
  \item{file}{A file name for the allele report. By default a sequential filename is created to avoid accidental overwriting.}
}
\description{Outputs a docx summary of the data inputs, to assist the user in choosing parameter inputs for a full evaluation}
\details{
The allele report summarises alleles present in both the crime scene profile 
and the reference profiles, prior to the main statistical evaluation. Parameter
inputs for nUnknowns and doDropin are suggested. 
The report is placed in the directory specified by admin$outputPath.
}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"hammer")

# File paths and case name for allele report
admin = pack.admin.input(
            cspFile = file.path(datapath, 'hammer-CSP.csv'),
            refFile = file.path(datapath, 'hammer-reference.csv'),
            caseName = "hammer",
	    kit= "SGMplus"
             )

# Next we generate an allele report
allele.report(admin)}
}
