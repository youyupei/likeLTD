## Test unit 'scenario'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_scenario.R"
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}

.tearDown <-
function () {
	## Specific actions for svUnit: clean up context
	if ("package:svUnit" %in% search()) {
		.Log$..Unit <- ""
		.Log$..File <- ""
		.Log$..Obj <- ""
		.Log$..Tag <- ""
		.Log$..Msg <- ""
		rm(..Test, envir = .Log)
	}
}


###############################################################
# Then data functions
###############################################################

ref.data.path <- function() {
  path = Reduce(file.path, c("extdata", "hammer", "hammer-reference.csv"))
  system.file(path, package="likeLTD")
}
csp.data.path <- function() {
  path = Reduce(file.path, c("extdata", "hammer", "hammer-CSP.csv"))
  system.file(path, package="likeLTD")
}


###################################
# Finally, the unit-test themselves
###################################
test_determine.dropout <- svTest(function() {

  cspProfile = read.csp.profile(csp.data.path())
  knownProfiles = read.known.profiles(ref.data.path())
  
  result = determine.dropout(knownProfiles, cspProfile)
  check = list("Suspect"=TRUE, "Victim 1"=TRUE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))

  # construct fake matrix with one replicate containing alleles from 1 and 2.
  other = mapply(union, knownProfiles[2, colnames(cspProfile), drop=F],
                 knownProfiles[1, colnames(cspProfile), drop=F])
  other = matrix(other, nrow=1)
  colnames(other) = colnames(cspProfile)
  result = determine.dropout(knownProfiles, other)
  check = list("Suspect"=FALSE, "Victim 1"=FALSE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))

  # construct fake matrix with one replicate containing alleles from 1 and 2
  # and another replicate containing alleles from 1 and 3.
  other2 = mapply(union, knownProfiles[3, colnames(cspProfile), drop=F],
                  knownProfiles[1, colnames(cspProfile), drop=F])
  other = rbind(other, other2)
  result = determine.dropout(knownProfiles, other)
  check = list("Suspect"=FALSE, "Victim 1"=TRUE, "Victim 2"=TRUE)
  checkTrue(result == check)
  checkTrue(names(result) == names(check))
})

test_missing.alleles <- svTest(function() {
  
  alleleDb   = list(D2=matrix(c(2, 2), nrow=14, ncol=2))
  row.names(alleleDb$D2) = c("10", "11", "12", "13", "14", "14.2", "15", "16",
                             "17", "18", "19", "20", "21", "22")
  cspProfile = matrix(list("D2"=list(c("10", "16"), c("14", "16", "17"))))
  noDropoutProfiles = matrix(list("D2"=list(c("10", "16"), c("14", "16"))))
  
  if(! "missing.alleles" %in% ls(.GlobalEnv))
    missing.alleles <- getFromNamespace("missing.alleles", "likeLTD")
  checkException(missing.alleles(alleleDb, cspProfile, noDropoutProfiles))

  colnames(cspProfile) = c("D2")
  colnames(noDropoutProfiles) = c("D2")
  result = missing.alleles(alleleDb, cspProfile, noDropoutProfiles)
  checkEquals(alleleDb, result)

  cspProfile = matrix(list("D2"=list(c("9", "16"), c("14", "16", "17"))))
  noDropoutProfiles = matrix(list("D2"=list(c("9", "16"), c("14", "16"))))
  colnames(cspProfile) = c("D2")
  colnames(noDropoutProfiles) = c("D2")
  checkEquals(alleleDb, missing.alleles(alleleDb, cspProfile, noDropoutProfiles))

  noDropoutProfiles = matrix(list("D2"=list(c("10", "16"), c("14", "16"))))
  colnames(noDropoutProfiles) = c("D2")
  result = missing.alleles(alleleDb, cspProfile, noDropoutProfiles)
  checkEquals(nrow(alleleDb$D2)+1, nrow(result$D2))
  checkTrue("9" %in% rownames(result$D2))
  checkEquals(result$D2["9", ], c(1, 0))
})
