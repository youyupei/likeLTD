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




