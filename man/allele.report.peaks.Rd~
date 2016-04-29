\name{allele.report.peaks} 
\alias{allele.report.peaks} 
\title{likeLTD::allele.report.peaks} 
\usage{allele.report.peaks(admin,file=NULL, figRes=300, dropinThresh=3)}
\arguments{
  \item{admin}{List containing administration data for peak heights inputs, 
as packed by pack.admin.input.peaks()}
  \item{file}{A file name for the allele report. By default a sequential filename 
is created to avoid accidental overwriting.}
\item{figRes}{Resolution of figures in the allele report. Defaults to 300 DPI.}
\item{dropinThresh}{If the estimated DNA contribution of Q is more than 
\code{dropinThresh} times the estimated DNA contribution of any unknown (estimated
by k-means) then it will be suggested that that unknown contributor could be 
modelled as dropin instead.}
}
\description{Outputs a docx summary of the data inputs, to assist the user in choosing 
parameter inputs for a full evaluation}
\details{
The allele report summarises the crime scene profile 
and the reference profiles, prior to the main statistical evaluation. Parameter
inputs for nUnknowns and doDropin are suggested. 
The report is placed in the directory specified by admin$outputPath.
}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"laboratory")

# File paths and case name for allele report
admin = pack.admin.input.peaks(
            peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
            refFile = file.path(datapath, 'laboratory-reference.csv'),
            caseName = "Laboratory",
            detectionThresh = 20
             )

# generate allele report
allele.report.peaks(admin)}
}
