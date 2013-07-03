\name{output.report}
\alias{output.report}
\title{likeLTD::output.report}
\usage{output.report(admin,prosecutionHypothesis,defenseHypothesis, prosecutionResults,defenseResults)}
\arguments{
  \item{admin}{Generated from pack.admin.input()}
  \item{prosecutionHypothesis}{The output from prosecution.hypothesis().}
  \item{defenseHypothesis}{The output from defense.hypothesis().}
  \item{prosecutionResults}{The output from do.call(optim, prosecutionParameters), or multiRun(prosecutionHypothesis, nrun, ...).}
  \item{defenseResults}{The output from do.call(optim, defenseParameters), or multiRun(defenseHypothesis, nrun, ...)$bestResult.}
}
\value{.pdf}
\description{Output results pdf forlikeLTD.}
\details{

}
\seealso{pdf, output.names}
