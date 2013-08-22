\name{output.names} 
\alias{output.names} 
\title{likeLTD::output.names} 
\usage{output.names(admin, prosecutionHypothesis)}
\arguments{
  \item{admin}{Arguments packed by pack.admin.input().}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}

}
\value{Two strings, pdf.name and RData.name.} 
\description{Automatic names for pdf and RData output.}
\details{
  For use with output.report(), to generate an automatic name for the 
  pdf and RData outputs.
}
\seealso{output.report}
