\name{output.report}
\alias{output.report}
\title{likeLTD::output.report}
\usage{output.report(admin,prosecutionHypothesis,defenceHypothesis, 
                      prosecutionResults,defenceResults,file)}
\arguments{
  \item{admin}{List containing administration data, as packed by pack.admin.input()}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenceHypothesis}{The output from defence.hypothesis().}
  \item{prosecutionResults}{The output from do.call(DEoptim, prosecutionParameters).}
  \item{defenceResults}{The output from do.call(DEoptim, defenceParameters).}
  \item{file}{A file name for the pdf. By default a filename is created from the amalgamation of parameters in admin.}
}
\description{Outputs a pdf of results.}
\details{
Creates a pdf containing all results, placed in the directory specified by admin$outputPath.
Results include both the pre-evaluation assessment of the alleles from both the reference profiles and the crime scene profile,
and the post-evaluation results of the likelihoods and weight of evidence.
}
\seealso{pdf, output.names}
