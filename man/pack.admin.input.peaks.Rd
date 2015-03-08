\name{pack.admin.input.peaks} 
\alias{pack.admin.input.peaks} 
\title{likeLTD::pack.admin.input.peaks} 
\usage{pack.admin.input.peaks(peaksFile, refFile, caseName='dummy',
                        databaseFile=NULL, kit=NULL, linkageFile=NULL,
                         outputPath=getwd())}
\arguments{
  \item{peaksFile}{Crime scene profile with peak height data}
  \item{refFile}{Reference profiles}
  \item{caseName}{Name of the current case}
  \item{databaseFile}{Path to the allele database}
  \item{kit}{Database to use if \code{databaseFile} not specified}
  \item{linkageFile}{Path to file with recombination rates for linked loci}
  \item{outputPath}{Path where the output should be stored. Defaults to current
  working directory.}
}
\description{Packs and verifies administrative information. Can be used as partial
input for \code{defence.hypothesis.peaks} and \code{prosecution.hypothesis.peaks}}




