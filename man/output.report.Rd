\name{output.report}
\alias{output.report}
\title{likeLTD::output.report}
\usage{output.report(admin,prosecutionHypothesis,defenceHypothesis, prosecutionResults,defenceResults)}
\arguments{
  \item{admin}{Generated from pack.admin.input()}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenceHypothesis}{The output from defence.hypothesis().}
  \item{prosecutionResults}{The output from do.call(optim, prosecutionParameters), or multiRun(prosecutionHypothesis, nrun, ...).}
  \item{defenceResults}{The output from do.call(optim, defenceParameters), or multiRun(defenceHypothesis, nrun, ...)$bestResult.}
}
\value{.pdf}
\description{Output results pdf forlikeLTD.}
\details{

}
\seealso{pdf, output.names}
