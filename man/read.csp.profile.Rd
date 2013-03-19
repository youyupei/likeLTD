\name{read.csp.profile} 
\alias{read.csp.profile} 
\title{likeLTD::read.csp.profile} 
\usage{read.csp.profile(path)}
\arguments{
  \item{path}{Path to the crime scene profile. }
}
\value{
  A matrix where columns are loci and rows are replicates. Each element is a
  vector of character, and each character is the name of an allele
  \emph{certainly} present in the CSP for that loci and replicate. For
  simplicity, loci for which all entries are NA are removed.
}
\description{Reads the Crime Scene Profile from a CSV file.}
\details{The input is a file in the CSV format (comma-separated values). 
It should have the following form:

  \tabular{llllccc}{
    Stain \tab Profiling system \tab  Plate/Run \tab Allelic/Uncertain \tab
      Locus \tab Locus \tab \dots \cr
    IGNORED \tab IGNORED \tab IGNORED \tab Allelic \tab "13,14,15" \tab "16,16" 
      \tab \dots\cr
    IGNORED \tab IGNORED \tab IGNORED \tab Uncertain \tab \tab "17.2"
      \tab \dots\cr
  }

Columns are separted by commas (not present in the table above). The first
line, containing the name of the columns need be present. 
The first three columns are ignored for the purpose of the \pkg{likeLTD}
package. They are present (and expected) merely for convenience. The fourth
column indicates whether a line contains certain or uncertain alleles.
"Uncertain" implies the latter, anything else implies the former. Beware of the
capitalization.  The next columns are for loci. "Locus" should be replaced with
the actual locus name in the first line. Other lines should contain the
comma-separated names of the alleles present in the CSP, within quotation
marks. Entries can be empty. The lines labeled "Uncertain" are ignored by this
function.
}
\seealso{read.unc.profile, read.known.profiles}
