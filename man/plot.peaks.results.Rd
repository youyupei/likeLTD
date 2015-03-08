\name{plot.peaks.results} 
\alias{plot.peaks.results} 
\title{likeLTD::plot.peaks.results} 
\usage{plot.peaks.results(hyp,res,replicate=1,fileName="peakHeights.pdf")}
\arguments{
    \item{hyp}{Hypothesis used to generate parameters for optimisation.}
    \item{res}{Either prosecution or defence results from \code{evaluate.peaks}
    e.g. \code{results$Pros} or \code{results$Def}.}
    \item{replicate}{Which replicate results to plot.}
    \item{fileName}{Output file name for plot.}
}
\value{
pdf file
} 
\description{
Plot distribution of fitted peaks along with actual peak heights. 
}
\details{
  CSP peak heights for a single replicate are plotted, with boxplots representing 
  the distribution for each peak height estimated through optimisation.
}
\seealso{plot.CSP.heights}




