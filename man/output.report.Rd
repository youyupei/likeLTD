\name{output.report}
\alias{output.report}
\title{likeLTD::output.report}
\usage{output.report(admin,prosecutionHypothesis,defenceHypothesis, prosecutionResults,defenceResults,file)}
\arguments{
  \item{admin}{Generated from pack.admin.input()}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenceHypothesis}{The output from defence.hypothesis().}
  \item{prosecutionResults}{The output from do.call(DEoptim, prosecutionParameters).}
  \item{defenceResults}{The output from do.call(DEoptim, defenceParameters).}
  \item{file}{A file name for the pdf. By default a filename is created from the amalgamation of paramters in admin.}
}
\value{}
\description{Outputs a pdf of results.}
\details{
  Creates a pdf containing all results, placed in the directory specified by admin$outputPath.
}
\seealso{pdf, output.names}
