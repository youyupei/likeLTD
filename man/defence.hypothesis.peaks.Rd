\name{defence.hypothesis.peaks}
\alias{defence.hypothesis.peaks}
\title{likeLTD::defence.hypothesis.peaks}
\description{Helper function to create the input for the defence using peak height data.}
\usage{defence.hypothesis.peaks(peaksFile, refFile, ethnic='NDU1', nUnknowns=0,
                        adj=1e0, fst=0.03, databaseFile=NULL, 
                        linkageFile=NULL,
                        detectionThresh=20, doDropin=FALSE, 
                        doDoubleStutter=TRUE,doOverStutter=TRUE,
                        combineRare=TRUE, rareThreshold=1,
                        kit=NULL, relationship=0, ...)
}
\arguments{
  \item{peaksFile}{Path to the crime scene profile.}
  \item{refFile}{Path to the known profiles.}
  \item{ethnic}{Ethnicity, e.g. subgroup within the allele frequency database.}
  \item{nUnknowns}{Number of unknown contributors for which to perform
                   calculation.}
  \item{adj}{Allele frequency adjustment parameter.}
  \item{fst}{Allele fraction adjustment \eqn{F_{ST}}.}
  \item{databaseFile}{Path to the allele database. If NULL, then defaults to
                      the NGMSelect database provided with \code{likeLTD}. Must include
                      longest uninterrupted sequence (LUS) values for alleles.}
  \item{linkageFile}{Path to recombination rate information. If NULL then defaults
                    to the linkage file provided with \code{likeLTD}. Only used
                    if \code{relationship} is not 0 or 1 i.e. when Q and X are 
                    closely related but not parent/offspring.}
  \item{relationship}{Specified relationship between Q and X. Can take values of 0=unrelated,
1=parent/offspring, 2=siblings, 3=uncle/nephew, 4=half-uncle/half-nephew, 5=cousins, 
6=grandparent/grandchild, 7=half-siblings.}
  \item{detectionThresh}{Detection threshold for peaks. Can be a single value, or a named list 
		containing one value per locus.}
  \item{doDropin}{Whether or not to model drop-in. Note dropin is not currently possible
                    with the peak heights model.}
  \item{doDoubleStutter}{Logical. Whether or not to model double stutter.}
  \item{doOverStutter}{Logical. Whether or not to model over stutter.}                  
  \item{combineRare}{Whether or not to combine rare unobserved alleles into a single
                    allele.}
  \item{rareThreshold}{If \code{combineRare=TRUE}, this is the probability threshold
                        below which an allele is classed as rare, and therefore 
                        combined with other rare unobserved alleles.}
  \item{kit}{Parameter specifying which allele database supplied with likeLTD to use if
            \code{linkageFile} is not specified. Possibilities are "DNA17", "Identifiler",
            "SGMplus", "NISTidentifiler" and "NISTglobalfiler".}
  \item{...}{Other parameters to determine how to perform calculations.}
}
\value{A list of named input parameters, or hypothesis, suitable for the
defence.}
\details{It loads the CSP, known profiles, and allele database from file. It
removes the queried individual from the known profiles. It increments the
number of unknown contributors by one (to make up for the queried individual).
} 
\seealso{prosecution.hypothesis.peaks}
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
}
}

