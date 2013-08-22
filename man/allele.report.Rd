\name{allele.report} 
\alias{allele.report} 
\title{likeLTD::allele.report} 
\usage{allele.report(admin,file)}
\arguments{
  \item{admin}{List containing administration data, as packed by pack.admin.input()}
  \item{file}{A file name for the pdf. By default a filename is created from the 
  amalgamation of parameters in admin.}
}
\description{Outputs a pdf summary of the data inputs, before a hypothesis is chosen}
\details{
The allele report summarises alleles present in both the crime scene profile 
and the reference profiles, prior to the main statistical evaluation. 
The report is placed in the directory specified by admin$outputPath.
}
