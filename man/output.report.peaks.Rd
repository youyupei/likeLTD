\name{output.report.peaks}
\alias{output.report.peaks}
\title{likeLTD::output.report.peaks}
\usage{output.report(prosecutionHypothesis,defenceHypothesis, 
                      results,file=NULL)}
\arguments{
  \item{prosecutionHypothesis}{The output from \code{prosecution.hypothesis.peaks}.}
  \item{defenceHypothesis}{The output from \code{defence.hypothesis.peaks}.}
  \item{results}{The output from \code{evaluate.peaks}.}
  \item{file}{A file name for the output report. By default a filename is 
created to avoid over writing a previous file.}
}
\description{Outputs a docx report of results.}
\details{
Creates a docx containing all results, placed in the directory specified 
by admin$outputPath.
Results include both the pre-evaluation assessment of the crime-stain and reference profiles, and the post-evaluation results of the likelihoods 
and weight of evidence.
}

