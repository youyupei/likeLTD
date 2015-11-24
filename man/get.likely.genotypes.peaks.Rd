\name{get.likely.genotypes.peaks}
\alias{get.likely.genotypes.peaks}
\title{likeLTD::get.likely.genotypes.peaks}
\usage{get.likely.genotypes.peaks(hypothesis,params,results,
posterior=FALSE,joint=FALSE,prob=ifelse(joint==FALSE,0.1,0.05))}
\arguments{
\item{hypothesis}{Hypothesis object created by either
\code{defence.hypothesis.peaks} or \code{prosecution.hypothesis.peaks}
that was used for optimisation with \code{evaluate.peaks}.}
\item{params}{Parameters object created by 
\code{optimisation.params.peaks} that was used for optimisation
with \code{evaluate}.}
\item{results}{Either prosecution or defence results returned by 
\code{evaluate.peaks} e.g. \code{results$Pros} or \code{results$Def}.}
\item{posterior}{Logical indicating whether to return all genotype
probabilities, rather than just the most likely.}
\item{joint}{Logical indicating whether or not to return joint
genotypes and probabilities. If FALSE, marginal genotypes and
probabilities are returned instead.}
\item{prob}{Probability threshold for single-locus genotype
probabilities. Defaults to 0.1 if returning marginal 
probabilities, and 0.05 if returning joint probabilities.}
}
\value{
\item{locusSpecific}{Locus genotypes and probabilities which are 
greater than \code{prob} for each contributor of \code{joint=FALSE}, 
or for all contributors if \code{joint=TRUE}.}
\item{topGenotype}{Most likely whole-profile genotype for all
contributors if \code{joint=TRUE} or for each contributor separately if
\code{joint=FALSE}. The probability of each genotype is also given.}
}
\description{Creates a list of the most likely genotypes at each locus,
and the most likely whole-profile genotype for peak height data.}
\details{Either joint or marginal genotypes and genotype probabilities 
are given. Locus-specific genotypes are only given if their probabilty
exceeds \code{prob}. The most likely whole-profile genotype is given,
regardless of the probability threshold at each locus. Joint 
probabilities give the probability of a multi-contributor genotype,
whereas marginal probabilities give the probability of a single
contributor, summing over all the possible genotypes for all other 
contributors.}
\seealso{defence.hypothesis.peaks, prosecution.hypothesis.peaks, 
optimisation.params.peaks,evaluate.peaks}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"laboratory")

# File paths and case name for allele report
admin = pack.admin.input.peaks(
            peaksFile = file.path(datapath, 'laboratory-CSP.csv'),
            refFile = file.path(datapath, 'laboratory-reference.csv'),
            caseName = "Laboratory",
            detectionThresh = 20
             )

# Enter arguments
args = list(
        nUnknowns = 1
        )

# Create hypotheses
hypP = do.call(prosecution.hypothesis.peaks, append(admin,args))
hypD = do.call(defence.hypothesis.peaks, append(admin,args))

# Get parameters for optimisation
paramsP = optimisation.params.peaks(hypP)
paramsD = optimisation.params.peaks(hypD)

# reduce number of iterations for demonstration purposes
paramsP$control$itermax=25
paramsD$control$itermax=25

# Run optimisation
# n.steps and converge set for demonstration purposes
results = evaluate.peaks(paramsP, paramsD, n.steps=1, 
      converge=FALSE)

# get most likely marginal genotypes under defence
get.likely.genotypes.peaks(hypD,paramsD,
      results$Def)

# get most likely joint genotypes under defence
gensJoint = get.likely.genotypes.peaks(hypD,paramsD,
      results$Def,joint=TRUE)

# get posterior likelihoods for all genotype combinations
  gensPosterior = get.likely.genotypes.peaks(hypD,paramsD,
      results$Def,posterior=TRUE)
}
}
