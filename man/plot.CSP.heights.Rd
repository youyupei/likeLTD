\name{plot.CSP.heights} 
\alias{plot.CSP.heights} 
\title{likeLTD::plot.CSP.heights} 
\usage{plot.CSP.heights(cspFile,refFile=NULL,detectThresh=NULL,
        stutterThresh=NULL,replicate=1)}
\arguments{
  \item{cspFile}{Path to crime scene profile file, with peak height data.}
  \item{refFile}{Path to file containing reference profiles. If \code{NULL} 
  peaks are not coloured according to the reference profiles.}
  \item{detectThresh}{Detection threshold for peaks. If \code{NULL} the 
  detection threshold is not plotted on the graph.}
  \item{stutterThresh}{Threshold for suggesting an peak is non-allelic. If
  \code{NULL} no information is presented on the graph suggesting allelic/non-allelic 
  calls.}
  \item{replicate}{Replicate to plot.}
}
\value{
Plot.
} 
\description{
Plot one replicate from a peak heights CSP.
}
\details{
  If \code{refFile} is specified, the plotted peaks are coloured according to which 
  reference individuals possess that allele. If \code{detectThresh} is specified, a 
  horizontal line is plotted indicating the threshold. If \code{stutterThresh} is 
  specified, the labels for each peak are coloured by a crude estimation of which
  peaks are allelic and which are not.
}
\seealso{allele.report.peaks}



