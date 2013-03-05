## Test unit 'genetics'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_genetics.R"
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

ref.data <- function() {
  # Reference profile used throughout the tests. 
  ref = list( D3   = c("14,16", "16,16", "15,17"),
              vWA  = c("15,19", "15,16", "16,19"),
              D16  = c("11,14", "13,13", "12,13"),
              D2   = c("24,25", "20,20", "18,25"),
              D8   = c("12,13", "11,15", "11,13"),
              D21  = c("28,31", "29,30", "29,30"),
              D18  = c("14,17", "17,17", "15,17"),
              D19  = c("15.2,17.2", "12,14", "14,14"),
              TH01 = c("9,9.3", "6,8", "6,7"),
              FGA  = c("22,23", "22,25", "20,22") )
  ref = data.frame( ref, row.names=c('Suspect', 'Victim 1', 'Victim 2'),
                    stringsAsFactors=FALSE )
  return(ref)
} 

ethnic.database.data <- function() {
  # Creates ethnic frequency data for tests.
  d3 = c( 1, 5, 115, 232, 216, 170, 123, 12, # first column
          -93.9413043478, -89.9413043478, -85.9413043478, -81.9413043478,
          -77.9413043478, -73.9413043478, -69.9413043478, -65.9413043478 )
  d3 = matrix(d3, , 2)
  row.names(d3) <- c("12", "13", "14", "15", "16", "17", "18", "19")
  d18 = c( 12, 5, 138, 113, 155, 1, 124, 107, 90, 69, 36, 15, 7, 2, # first column
           70.0586956522, 74.0586956522, 78.0586956522, 82.0586956522,
           86.0586956522, 88.0586956522, 90.0586956522, 94.0586956522,
           98.0586956522, 102.0586956522, 106.0586956522, 110.0586956522,
           114.0586956522, 118.0587 )
  d18 = matrix(d18, , 2)
  row.names(d18) <- c("10", "11", "12", "13", "14", "14.2", "15", "16", "17",
                      "18", "19", "20", "21", "22")
  d19 = c( 76, 194, 11, 334, 13, 155, 33, 36, 15, 4, 2, 1,
           -90.9413043478, -86.9413043478, -84.9413043478, -82.9413043478, -80.9413043478,
           -78.9413043478, -76.9413043478, -74.9413043478, -72.9413043478, -70.9413043478,
           -64.9413043478, -60.9413043478 ) 
  d19 = matrix(d19, , 2)
  row.names(d19) <- c("12", "13", "13.2", "14", "14.2", "15", "15.2", "16",
                      "16.2", "17", "18.2", "19.2")

  return(list(D3=d3, D18=d18, D19=d19))
}

internal.representation.data = function() {
  # Internal represenation data as reconstructed by hand.
  # 
  # Used for testing both the construction of the internal representation and
  # the pack.genetic.input test.
  result = list( D3=  list( list(csp=c("14", "16"), unc=""),
                            list(csp=c("14", "16"), unc="") ),
                 vWA= list( list(csp=c("15", "16", "19"), unc=""),
                            list(csp=c("15", "16", "17", "19"), unc="") ),
                 D16= list( list(csp=c("11", "13", "14"), unc=""),
                            list(csp=c("11", "13", "14"), unc="") ),
                 D2=  list( list(csp=c("20", "23", "24", "25"), unc=""),
                            list(csp=c("20", "24", "25"), unc="") ),
                 D8=  list( list(csp=c("11", "12", "13", "15"), unc=""),
                            list(csp=c("11", "12", "13", "15"), unc="") ),
                 D21= list( list(csp=c("28", "31"), unc=""),
                            list(csp=c("28", "29", "30", "31", "31.2"), unc="") ),
                 D18= list( list(csp=c(""), unc=""),
                            list(csp=c("13", "14", "16", "17"), unc="") ),
                 D19= list( list(csp=c("12", "13", "15.2", "17.2"), unc=""),
                            list(csp=c("12", "13", "14", "15.2", "17.2"), unc="") ),
                 TH01=list( list(csp=c("6", "8", "9", "9.3"), unc=""),
                            list(csp=c("6", "8", "9", "9.3"), unc="") ),
                 FGA= list( list(csp=c("22"), unc=""), 
                            list(csp=c("22", "23", "25"), unc="") ) )
  return(result)
}

known.alleles.data <- function() {
  # Result from known allele call. 
  check <- list( D3   = c("14", "16", "16", "16", "15", "17"), 
                 vWA  = c("15", "19", "15", "16", "16", "19"), 
                 D16  = c("11", "14", "13", "13", "12", "13"), 
                 D2   = c("24", "25", "20", "20", "18", "25"), 
                 D8   = c("12", "13", "11", "15", "11", "13"), 
                 D21  = c("28", "31", "29", "30", "29", "30"), 
                 D18  = c("14", "17", "17", "17", "15", "17"), 
                 D19  = c("15.2", "17.2", "12", "14", "14", "14"), 
                 TH01 = c("9", "9.3", "6", "8", "6", "7"), 
                 FGA  = c("22", "23", "22", "25", "20", "22") )
  return(data.frame(check, stringsAsFactors=FALSE))
}

###################################
# Finally, the unit-test themselves
###################################
test_ethnic.database <- svTest(function() {
  # Test ethnic frequency function
  if(! "ethnic.database" %in% ls(.GlobalEnv))
    ethnic.database = getFromNamespace("ethnic.database", "likeLTD")
  result = ethnic.database('EA1')

  check <- ethnic.database.data()
  checkEquals(result$D3, check$D3)
  checkEquals(result$D18, check$D18)
  checkEquals(result$D19, check$D19)
})


test_all.genotypes.per.locus <- svTest(function() {
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus = getFromNamespace("all.genotypes.per.locus", "likeLTD")
  result = all.genotypes.per.locus(2)
  checkEquals(result, t(matrix(c(1, 1, 2, 1, 2, 2), nrow=3)))
  
  result = all.genotypes.per.locus(2, 2)
  check = matrix( c(1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1,
                    2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2), nrow=9 )
  checkEquals(result, t(check))
  
  result = all.genotypes.per.locus(2, 3)
  check = matrix( c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
                    2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2,
                    2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1,
                    1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2,
                    1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 2,
                    2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1,
                    2, 2), ncol=6 )
  checkEquals(result, t(check))

  result = all.genotypes.per.locus(3)
  check = matrix(c(1, 1, 1, 2, 2, 3, 1, 2, 3, 2, 3, 3), ncol=2,)
  checkEquals(result, t(check))
  
  result = all.genotypes.per.locus(3, 2)
  check = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2,
                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1,
                   1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,
                   3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 3, 1, 1,
                   1, 2, 2, 3, 1, 1, 1, 2, 2, 3, 1, 1, 1, 2, 2, 3, 1, 1, 1, 2,
                   2, 3, 1, 1, 1, 2, 2, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3, 2, 3, 3,
                   1, 2, 3, 2, 3, 3, 1, 2, 3, 2, 3, 3, 1, 2, 3, 2, 3, 3, 1, 2,
                   3, 2, 3, 3), ncol=4) 
  checkEquals(result, t(check))
})

test_possible.genotypes <- svTest(function() {
  alleleNames = c("one", "two", "three", "four", "five")
  cspPresence = matrix(c(TRUE, FALSE, FALSE, TRUE, FALSE), ncol=1)
  profPresence = matrix(c(TRUE, TRUE, TRUE, FALSE, FALSE), ncol=1)
  missingReps = FALSE
  # Basic trial
  if(! "possible.genotypes" %in% ls(.GlobalEnv))
    possible.genotypes <- getFromNamespace("possible.genotypes", "likeLTD")
  result <- possible.genotypes(cspPresence, profPresence,
                               missingReps, alleleNames, 1, FALSE)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 5)
  for(i in 1:nrow(result)) checkTrue(4 %in% result[i,])
  checkTrue(all(result[, 1] <= result[, 2]))

  # Try with csp matrix.
  cspPresence = t(matrix(c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           TRUE, FALSE, FALSE), nrow=2))
  result <- possible.genotypes(cspPresence, profPresence,
                               missingReps, alleleNames, 1, FALSE)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 5)
  for(i in 1:ncol(result)) checkTrue(4 %in% result[, i])
  checkTrue(all(result[1, ] <= result[2, ]))

  # Try with too large CSP vs contribs
  cspPresence = matrix(c(TRUE, TRUE, TRUE, TRUE, TRUE), ncol=1)
  profPresence = matrix(c(TRUE, FALSE, FALSE, FALSE, FALSE), ncol=1)
  checkException( possible.genotypes(cspPresence, profPresence,
                                     missingReps, alleleNames, 1, FALSE) )
  # Try with exact match
  result <- possible.genotypes(cspPresence, profPresence,
                              missingReps, alleleNames, 2, FALSE)
  checkEquals(result, matrix(0, 1, 0))

  # Try with two contributors 
  profPresence = matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE), ncol=1)
  result <- possible.genotypes(cspPresence, profPresence,
                              missingReps, alleleNames, 2, FALSE)
  checkTrue(ncol(result) == 24)
  for(i in 1:ncol(result)) checkTrue(all(3:5 %in% result[, i]))

  # Try with non-matrix input.
  cspPresence = c(TRUE, TRUE, TRUE, TRUE, TRUE)
  checkException( possible.genotypes(cspPresence, profPresence,
                                     missingReps, alleleNames, 1, FALSE) )
  checkException( possible.genotypes(cspPresence, profPresence,
                                     missingReps, alleleNames, 1, TRUE) )
  checkException( possible.genotypes(profPresence, cspPresence,
                                     missingReps, alleleNames, 1, FALSE) )
  checkException( possible.genotypes(profPresence, cspPresence,
                                     missingReps, alleleNames, 1, TRUE) )
})

test_possible.profiles.with.dropin <- svTest(function() {
  alleleNames = c("one", "two", "three", "four", "five")
  cspPresence = matrix(ncol=0, nrow=length(alleleNames))
  profPresence = matrix(ncol=0, nrow=length(alleleNames))
  missingReps = c()

  # Basic trial
  if(! "possible.genotypes" %in% ls(.GlobalEnv))
    possible.genotypes <- getFromNamespace("possible.genotypes", "likeLTD")
  result <- possible.genotypes(cspPresence, profPresence,
                               missingReps, alleleNames, 1, TRUE)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 15)
  checkTrue(nrow(result) == 2)
  checkEquals(result, t(combinations(5, 2, rep=TRUE)))

  # Try with two possible.profiles.
  result <- possible.genotypes(cspPresence, profPresence,
                              missingReps, alleleNames, 2, TRUE)
  checkTrue(ncol(result) == (5*3)^2)
  checkTrue(nrow(result) == 4)
  checkEquals(t(result), unique(t(result)))

  # Try with three possible.profiles.
  result <- possible.genotypes(cspPresence, profPresence,
                              missingReps, alleleNames, 3, TRUE)
  checkTrue(ncol(result) == (5*3)^3)
  checkTrue(nrow(result) == 6)
  checkEquals(t(result), unique(t(result)))
})

test_known.epg.per.locus = svTest(function() {

  profiles = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0), ncol=3) > 0
  # NA is to make sure we can pass vectors that are too long, e.g. with
  # unprofiled contributors.
  degradation = c(0.001, 0.002, 0.003, NA) 
  relContrib  = c(0.01, 0.02, 0.03, NA)
  database = matrix(c(1.316143e-02, 5.483929e-03, 1.513564e-01,
                      1.239368e-01, 1.907064e-01, 1.096786e-03,
                      1.360014e-01, 1.173561e-01, 1.194153e-01,
                      7.567822e-02, 3.948429e-02, 1.645179e-02,
                      7.677500e-03, 2.193571e-03, 7.005870e+01,
                      7.405870e+01, 7.805870e+01, 8.205870e+01,
                      8.605870e+01, 8.805870e+01, 9.005870e+01,
                      9.405870e+01, 9.805870e+01, 1.020587e+02,
                      1.060587e+02, 1.100587e+02, 1.140587e+02,
                      1.180587e+02), ncol=2)

  if(! "known.epg.per.locus" %in% ls(.GlobalEnv))
    known.epg.per.locus <- getFromNamespace("known.epg.per.locus", "likeLTD")
  result <- known.epg.per.locus(relContrib, degradation, database[, 2],
                                profiles)
  checkTrue(all(result[!rowSums(profiles)] == 0))
  checkEquals( result[rowSums(profiles) > 0][[1]],
               relContrib[1] * (1.0 + degradation[1])^-database[5, 2])
  checkEquals( result[rowSums(profiles) > 0][[2]], 
               relContrib[3] * (1.0 + degradation[3])^-database[7, 2] )
  checkEquals( result[rowSums(profiles) > 0][[3]], 
               relContrib[1] * (1.0 + degradation[1])^-database[9, 2] 
               + 2.0 * relContrib[2] * (1.0 + degradation[2])^-database[9, 2] 
               + relContrib[3] * (1.0 + degradation[3])^-database[9, 2]  )
})

test_selective.col.prod.mat.mat= svTest(function() {

  if(! "selective.col.prod" %in% ls(.GlobalEnv))
    selective.col.prod <- getFromNamespace("selective.col.prod", "likeLTD")
  
  condition = matrix(FALSE, ncol=6, nrow=3)
  input     = matrix(2, ncol=6, nrow=3)

  # Checks all False
  checkEquals(selective.col.prod(condition, input), 1)

  # Checks full rows can be unselected.
  condition = matrix(c(TRUE, FALSE), ncol=6, nrow=3, byrow=TRUE)
  checkEquals(selective.col.prod(condition, input), rep(c(2^3, 1), 3))

  # Checks more random selection.
  condition = matrix(c(TRUE, FALSE), ncol=6, nrow=3)
  checkEquals(selective.col.prod(condition, input), rep(c(4, 2), 3))

  # Check that this is not a matrix-vector selection
  input[1:2, 5:6] = 3 
  input[3, 5:6] = 4 
  input[3, 1:4] = 5
  checkEquals(selective.col.prod(condition, input), c(10, 2, 10, 2, 12, 3))

  condition = matrix(TRUE, ncol=7, nrow=3)
  checkException(selective.col.prod(condition, input))
  condition = matrix(TRUE, ncol=6, nrow=4)
  checkException(selective.col.prod(condition, input))
})

test_selective.col.prod.vec.mat= svTest(function() {

  if(! "selective.col.prod" %in% ls(.GlobalEnv))
    selective.col.prod <- getFromNamespace("selective.col.prod", "likeLTD")
  
  condition = rep(FALSE, 3)
  input     = matrix(2, ncol=6, nrow=3)

  # Checks all False
  checkEquals(selective.col.prod(condition, input), 1)

  condition = rep(TRUE, 3)
  checkEquals(selective.col.prod(condition, input), rep(8, 6))

  condition = c(TRUE, FALSE, FALSE)
  checkEquals(selective.col.prod(condition, input), rep(2, 6))

  condition = c(TRUE, FALSE, TRUE)
  input[2, ] = 6
  checkEquals(selective.col.prod(condition, input), rep(4, 6))
  
  input[3, ] = 6
  checkEquals(selective.col.prod(condition, input), rep(12, 6))

  condition = c(TRUE, FALSE, TRUE, TRUE)
  checkException(selective.col.prod(condition, input))
})

test_selective.col.prod.mat.vec = svTest(function() {

  if(! "selective.col.prod" %in% ls(.GlobalEnv))
    selective.col.prod <- getFromNamespace("selective.col.prod", "likeLTD")
  
  condition = matrix(FALSE, ncol=6, nrow=3)
  input     = rep(2, 3)

  # Checks all False
  checkEquals(selective.col.prod(condition, input), 1)

  # Checks full rows can be unselected.
  condition = matrix(c(TRUE, FALSE), ncol=6, nrow=3, byrow=TRUE)
  checkEquals(selective.col.prod(condition, input), rep(c(2^3, 1), 3))

  # Checks more random selection.
  condition = matrix(c(TRUE, FALSE), ncol=6, nrow=3)
  checkEquals(selective.col.prod(condition, input), rep(c(4, 2), 3))

  input = 1:3
  condition = matrix(c(TRUE, FALSE), ncol=6, nrow=3)
  checkEquals(selective.col.prod(condition, input), rep(c(3, 2), 3))

  checkException(selective.col.prod(condition, 1:4))
})

test_selective.col.prod.vec.vec = svTest(function() {

  if(! "selective.col.prod" %in% ls(.GlobalEnv))
    selective.col.prod <- getFromNamespace("selective.col.prod", "likeLTD")

  checkEquals(selective.col.prod(rep(FALSE, 3), rep(2, 3)), 1)
  checkEquals(selective.col.prod(rep(TRUE, 3), rep(2, 3)), rep(2, 3))
  checkEquals(selective.col.prod(c(TRUE, FALSE, TRUE), rep(2, 3)), c(2, 1, 2))

  checkException(selective.col.prod(rep(TRUE, 3), 1:4))
})

