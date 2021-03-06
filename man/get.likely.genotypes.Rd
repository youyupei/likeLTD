\name{get.likely.genotypes}
\alias{get.likely.genotypes}
\title{likeLTD::get.likely.genotypes}
\usage{get.likely.genotypes(hypothesis,params,results,
posterior=FALSE,joint=FALSE,prob=ifelse(joint==FALSE,0.1,0.05))}
\arguments{
\item{hypothesis}{Hypothesis object created by either
\code{defence.hypothesis} or \code{prosecution.hypothesis}
that was used for optimisation with \code{DEoptimLoop}.}
\item{params}{Parameters object created by 
\code{optimisation.params} that was used for optimisation
with \code{DEoptimLoop}.}
\item{results}{Results object returned by \code{DEoptimLoop}}
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
and the most likely whole-profile genotype.}
\details{Either joint or marginal genotypes and genotype probabilities 
are given. Locus-specific genotypes are only given if their probabilty
exceeds \code{prob}. The most likely whole-profile genotype is given,
regardless of the probability threshold at each locus. Joint 
probabilities give the probability of a multi-contributor genotype,
whereas marginal probabilities give the probability of a single
contributor, summing over all the possible genotypes for all other 
contributors.}
\seealso{defence.hypothesis, prosecution.hypothesis, 
optimisation.params,DEoptimLoop}
\examples{
\dontrun{
# datapath to example files
datapath = file.path(system.file("extdata", package="likeLTD"),"hammer")

# File paths and case name for allele report
admin = pack.admin.input(
            cspFile = file.path(datapath, 'hammer-CSP.csv'),
            refFile = file.path(datapath, 'hammer-reference.csv'),
            caseName = "hammer",
	    kit= "SGMplus"
             )

# Enter arguments
args = list(
        nUnknowns = 1,
        doDropin = FALSE,
        ethnic = "EA1",
        adj = 1,
        fst = 0.02,
        relatedness = c(0,0)
        )

# Create hypotheses
hypP = do.call(prosecution.hypothesis, append(admin,args))
hypD = do.call(defence.hypothesis, append(admin,args))

# Get parameters for optimisation
paramsP = optimisation.params(hypP)
paramsD = optimisation.params(hypD)

# Run optimisation
# n.steps set for demonstration purposes
results = evaluate(paramsP, paramsD, n.steps=1)

# get most likely marginal genotypes under defence
get.likely.genotypes(hypD,paramsD,results$Def)

# get most likely joint genotypes under defence
get.likely.genotypes(hypD,paramsD,results$Def,joint=TRUE)

# get full posterior likelihoods
get.likely.genotypes(hypD,paramsD,results$Def,posterior=TRUE)
}
}

