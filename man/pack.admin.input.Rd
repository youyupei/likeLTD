\name{pack.admin.input} 
\alias{pack.admin.input} 
\title{likeLTD::pack.admin.input} 
\usage{pack.admin.input(mixedFile, refFile, caseName='dummy',
                        databaseFile=NULL, outputPath=NULL,
                        checkFiles=TRUE)}
\arguments{
  \item{databaseFile}{File holding the allele afreq}
  \item{mixedFile}{Mixed crime scene profile}
  \item{refFile}{Reference profiles}
  \item{casename}{Name of the current case}
  \item{outputPath}{Path where the output should be stored}
  \item{checkFiles}{Whether to check file existence and such}
}
\description{Packs and verifies administrative information. Only used in allele
report.}
