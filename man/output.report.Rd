\name{output.report}
\alias{output.report}
\title{likeLTD::output.report}
\usage{output.report(prosecutionHypothesis=NULL,defenceHypothesis=NULL, 
                      prosecutionResults=NULL,defenceResults=NULL,file=NULL)}
\arguments{
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenceHypothesis}{The output from defence.hypothesis().}
  \item{prosecutionResults}{The output from do.call(DEoptim, prosecutionParameters).}
  \item{defenceResults}{The output from do.call(DEoptim, defenceParameters).}
  \item{file}{A file name for the output report. By default a filename is created to avoid over writing a previous file.}
}
\description{Outputs a docx report of results.}
\details{
Creates a docx containing all results, placed in the directory specified 
by admin$outputPath.
Results include both the pre-evaluation assessment of the alleles from both the reference 
profiles and the crime scene profile, and the post-evaluation results of the likelihoods 
and weight of evidence.
}
\seealso{R2DOC, R2DOCX}
