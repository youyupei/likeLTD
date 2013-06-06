\name{output.report} 
\alias{output.report} 
\title{likeLTD::output.report} 
\usage{output.report(admin,prosecutionHypothesis,defenseHypothesis, prosecutionResults,defenseResults,args,stain)}
\arguments{
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenseHypothesis}{The output from defense.hypothesis().}
  \item{prosecutionResults}{The output from do.call(optim, prosecutionParameters), or multiRun(prosecutionHypothesis, nrun, ...)$bestResult.}
  \item{defenseResults}{The output from do.call(optim, defenseParameters), or multiRun(defenseHypothesis, nrun, ...)$bestResult.}
  \item{args}{Arguments specified by user that were passed to defense.hypothesis() or prosecution.hypothesis().}
  \item{caseName}{String containing the name of the case.}
  \item{stain}{String containing the name of the stain.}
}
\value{Plot.} 
\description{Output results plots forlikeLTD.}
\details{
  Best used with pdf(width=8.25,height=11.5,file). Can use outputNames()$pdf.name to generate automatic name.
}
\seealso{pdf, output.names}
