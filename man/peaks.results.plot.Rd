\name{peaks.results.plot} 
\alias{peaks.results.plot} 
\title{likeLTD::peaks.results.plot} 
\usage{peaks.results.plot(hyp,res,replicate=1,toplot=NULL,fileName=NULL,...)}
\arguments{
    \item{hyp}{Hypothesis used to generate parameters for optimisation.}
    \item{res}{Either prosecution or defence results from \code{evaluate.peaks}
    e.g. \code{results$Pros} or \code{results$Def}.}
    \item{replicate}{Which replicate results to plot.}
    \item{toplot}{Integer vector indicating which loci to plot. If \code{NULL},
	all loci are plotted.}
    \item{fileName}{Output file name for plot. If \code{NULL} the plot is output to 
	a plot window instead.}
    \item{...}{Extra parameters to pass to \code{boxplot}.}
}
\value{
Pdf file or plot window.
} 
\description{
Plot distribution of fitted peaks along with actual peak heights. 
}
\details{
  CSP peak heights for a single replicate are plotted, with boxplots representing 
  the distribution for each peak height estimated through optimisation.
}
\seealso{plot.CSP.heights}
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

# plot fitted results under prosecution
peaks.results.plot(hypP,results$Pros)

# plot fitted results under defence
peaks.results.plot(hypD,results$Def)

# plot fitted results for first four loci
peaks.results.plot(hypP,results$Pros,toplot=1:4)
}
}



