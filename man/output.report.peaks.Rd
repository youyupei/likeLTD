\name{output.report.peaks}
\alias{output.report.peaks}
\title{likeLTD::output.report.peaks}
\usage{output.report.peaks(prosecutionHypothesis,defenceHypothesis, 
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

# generate output report
output.report.peaks(hypP,hypD,results)
}
}
