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
  \item{file}{A file name for the .pdf}
}
\value{}
\description{Output results pdf forlikeLTD.}
\details{
  \item{output}{Creates a pdf containing all results, placed in the directory specified by admin$outputPath.}
  \item{file}{If the file name for the pdf is omitted, by default a filename is created from the amalgamation of paramters in admin.}
}
\seealso{pdf, output.names}
