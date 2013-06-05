\name{outputNames} 
\alias{outputNames} 
\title{likeLTD::outputNames} 
\usage{outputNames(stain, args, prosecutionHypothesis, caseFolder)}
\arguments{
  \item{stain}{String containing the name of the stain.}
  \item{args}{Arguments specified by user that were passed to defense.hypothesis() or prosecution.hypothesis().}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{caseFolder}{String containing the path to the directory you want to output to.}
}
\value{Two strings, pdf.name and RData.name.} 
\description{Automatic names for pdf and RData output.}
\details{
  For use with outputReport(), to generate an automatic name for the pdf and RData outputs.
}
\seealso{outputReport}
