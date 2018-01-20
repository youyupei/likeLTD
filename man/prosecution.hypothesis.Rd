\name{prosecution.hypothesis}
\alias{prosecution.hypothesis}
\title{likeLTD::prosecution.hypothesis}
\usage{prosecution.hypothesis(cspFile, refFile, ethnic='NDU1', nUnknowns=0,
                        adj=1e0, fst=0.02, databaseFile=NULL, 
                        linkageFile=NULL,
                        doDropin=FALSE, 
                        combineRare=TRUE, rareThreshold=0.05,
                        kit=NULL, relationship=0, ...)
}
\arguments{
   \item{cspFile}{Path to the crime scene profile.}
  \item{refFile}{Path to the known profiles.}
  \item{ethnic}{Ethnicity, e.g. subgroup within the allele frequency database.}
  \item{nUnknowns}{Number of unknown contributors for which to perform
                   calculation.}
  \item{adj}{Allele frequency adjustment parameter.}
  \item{fst}{Allele fraction adjustment \eqn{F_{ST}}.}
  \item{databaseFile}{Path to the allele database. If NULL, then defaults to
                      the NGMSelect database provided with \code{likeLTD}.}
  \item{linkageFile}{Path to recombination rate information. If NULL then defaults
                    to the linkage file provided with \code{likeLTD}. Only used
                    if \code{relatedness}=\code{c(0.5,0.25)} i.e. when Q and X are 
                    siblings.}
  \item{relationship}{Specified relationship between Q and X. Can take values of 0=unrelated,
1=parent/offspring, 2=siblings, 3=uncle/nephew, 4=half-uncle/half-nephew, 5=cousins, 
6=grandparent/grandchild, 7=half-siblings.}
  \item{doDropin}{Whether or not to model drop-in.}
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
prosecution.}
\description{Helper function to create the input for the prosecution.}
\details{It loads the CSP, known profiles, and allele database from file. It
makes sure the relatedness factor is zero. It sets the reference individual to
the first queried individual, unless it is not subject to dropout. In that
case, the reference individual is the first individual with dropout. } 
\seealso{defence.hypothesis}
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
}
}
