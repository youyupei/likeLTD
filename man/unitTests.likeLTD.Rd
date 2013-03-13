\name{unitTests.likeLTD}
\alias{unitTests.likeLTD}
\title{ Unit tests for the package likeLTD }
\description{ Performs unit tests defined in this
package by running \code{example(unitTests.likeLTD)}.
Tests are in \code{runit*.R} files Located in the
'/unitTests' subdirectory or one of its
subdirectories ('/inst/unitTests' and subdirectories
in package sources).
}
\author{John Doe
(\email{john.doe@doelabs.org})}
\examples{
if (require(svUnit)) {
clearLog()
runTest(svSuite("package:likeLTD"))
## Possibly run other tests here...
errorLog()
}
}
\keyword{utilities}
