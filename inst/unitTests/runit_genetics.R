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
# Then two functions to help design unit tests + data functions
###############################################################

temporary.directory <- function(expr) {
  # Creates a temporary working directory.
  #
  # The temp directory is created. The working directory is set to it. The
  # expression is evaluated. Whatever happens, the directory should be
  # removed when we exit this function.
  startdir = getwd()
  tryCatch( { 
              directory = tempfile()
              dir.create(directory, recursive=TRUE)
              setwd(directory)
              ;expr;
            },
            finally = {
              try(unlink(directory, TRUE, TRUE), silent=TRUE)
              setwd(startdir)
            } )
}

checkNoException <- function(expr, msg="") {
  foundError <- FALSE
  tryCatch(expr, error=function(n) foundError <- TRUE)
  checkTrue(!foundError, msg=msg)
}

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

ethnic.frequencies.data <- function() {
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
test_known.alleles <- svTest(function()  {
  # Tests known alleles
 
  if(! "known.alleles" %in% ls(.GlobalEnv))
    known.alleles = getFromNamespace("known.alleles", "likeLTD")
  result <- known.alleles(c("Suspect", "Victim 1"), ref.data()) 
  check = list( D3   = c("14", "16", "16", "16"), 
                vWA  = c("15", "19", "15", "16"), 
                D16  = c("11", "14", "13", "13"), 
                D2   = c("24", "25", "20", "20"), 
                D8   = c("12", "13", "11", "15"), 
                D21  = c("28", "31", "29", "30"), 
                D18  = c("14", "17", "17", "17"), 
                D19  = c("15.2","17.2", "12", "14"), 
                TH01 = c("9", "9.3", "6", "8"), 
                FGA  = c("22", "23", "22", "25") )
  check = data.frame(check, stringsAsFactors=FALSE)
  checkEquals(result, check)

  result <- known.alleles(c("Suspect", "Victim 1", "Victim 2"), ref.data())
  check <- known.alleles.data()
  checkEquals(result, check)
})

test_check.dropouts <- svTest(function() {
  # Checks that we can figure out dropouts.
  cprofs = internal.representation.data()

  if(! "has.dropouts" %in% ls(.GlobalEnv))
    has.dropouts = getFromNamespace("has.dropouts", "likeLTD")
  checkTrue(has.dropouts("Suspect", ref.data(), cprofs))

  cprofs$FGA[[1]]$csp <- c("22", "23")
  checkTrue(has.dropouts("Suspect", ref.data(), cprofs))

  cprofs$D18[[1]]$csp <- NA
  checkTrue(!has.dropouts("Suspect", ref.data(), cprofs))

  cprofs$D18[[1]]$csp <- c("14", "17")
  checkTrue(!has.dropouts("Suspect", ref.data(), cprofs))
})

test_ethnic.frequencies <- svTest(function() {
  # Test ethnic frequency function
  if(! "ethnic.frequencies" %in% ls(.GlobalEnv))
    ethnic.frequencies = getFromNamespace("ethnic.frequencies", "likeLTD")
  result = ethnic.frequencies('EA1')

  check <- ethnic.frequencies.data()
  checkEquals(result$D3, check$D3)
  checkEquals(result$D18, check$D18)
  checkEquals(result$D19, check$D19)
})


test_add.missing.alleles <- svTest(function() {
  # Add test for missing alleles.
  frq <- ethnic.frequencies.data()
  cprofs <- internal.representation.data()[c("D3", "D18", "D19")]
  known <- known.alleles.data()

  if(! "add.missing.alleles" %in% ls(.GlobalEnv))
    add.missing.alleles = getFromNamespace("add.missing.alleles", "likeLTD")

  # Checks add from queried
  result = add.missing.alleles(frq, cprofs, known[1:2, ])
  unchanged = setdiff(names(result), c("D19"))
  checknrow <- function(n) nrow(result[[n]]) == nrow(frq[[n]])
  checkTrue(all(sapply(unchanged, checknrow)))
  checkEquals(nrow(result$D19), nrow(frq$D19) + 1)
  checkEquals(row.names(result$D19)[nrow(result$D19)], "17.2")
  checkEquals(result$D19[nrow(result$D19),], c(1, 0))

  # Checks add from csp 
  result = add.missing.alleles(frq, cprofs, NULL)
  checknrow <- function(n) nrow(result[[n]]) == nrow(frq[[n]])
  checkTrue(all(sapply(unchanged, checknrow)))
  checkEquals(nrow(result$D19), nrow(frq$D19) + 1)
  checkEquals(row.names(result$D19)[nrow(result$D19)], "17.2")
  checkEquals(result$D19[nrow(result$D19),], c(1, 0))

  # Checks no add from csp but in profiled
  result = add.missing.alleles(frq, cprofs, NULL, known[1:2, ])
  checknrow <- function(n) nrow(result[[n]]) == nrow(frq[[n]])
  checkTrue(all(sapply(names(result), checknrow)))
})

test_presence.matrices <- svTest(function() {
  # Test construction of presence matrices.
  frq <- ethnic.frequencies.data()
  frq$D19 = rbind(frq$D19, c(1, 0))
  row.names(frq$D19)[nrow(frq$D19)] <- "17.2"
  cprofs <- internal.representation.data()[c("D3", "D18", "D19")]
  d3  = t(matrix(c(0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0), ncol=2))
  d18 = t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                   1, 1, 0, 0, 0, 0, 0), ncol=2))
  d19 = t(matrix(c(1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), ncol=2))
  check = list(D3=d3, D18=d18, D19=d19)

  if(! "presence.matrices" %in% ls(.GlobalEnv))
    presence.matrices = getFromNamespace("presence.matrices", "likeLTD")
  matrices <- presence.matrices(frq, cprofs)
  checkEquals(matrices, check)

  # Now check adding an extra matrix at the end.
  knownAlleles <- known.alleles.data()[, c("D3", "D18", "D19")]
  d3  = t(matrix(c(0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0), ncol=2))
  d18 = t(matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
                   1, 2, 0, 0, 0, 0, 0), ncol=2))
  d19 = t(matrix(c(2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 1, 0, 2, 0, 0, 1, 0,
                   0, 0, 0, 0, 1), ncol=2))
  check = list(D3=d3, D18=d18, D19=d19)
  matrices <- presence.matrices(frq, cprofs, knownAlleles[3:4,])
  checkEquals(matrices, check)
})

test_adjust.frequencies <- svTest(function(){ 
  # Test frequency adjustments.
  frq <- ethnic.frequencies.data()[c("D3", "D18", "D19")]
  frq$D19 = rbind(frq$D19, c(1, 0))
  row.names(frq$D19)[nrow(frq$D19)] <- "17.2"
  knownAlleles <- known.alleles.data()[1:2, c("D3", "D18", "D19")]

  
  d3 = list(c(0.00109678574626, -93.9413043478),
            c(0.00548392873131, -89.9413043478),
            c(0.14683498970364, -85.9413043478),
            c(0.25445429313278, -81.9413043478),
            c(0.25761035007610, -77.9413043478),
            c(0.18645357686454, -73.9413043478),
            c(0.13490464679022, -69.9413043478),
            c(0.01316142895514, -65.9413043478))
  d3 <- t(data.frame(d3))
  row.names(d3) <- c("12", "13", "14", "15", "16", "17", "18", "19")

  if(! "adjust.frequencies" %in% ls(.GlobalEnv))
    presence.matrices = getFromNamespace("adjust.frequencies", "likeLTD")
  result = adjust.frequencies(frq, knownAlleles)
  checkEquals(result$D3, d3)

  d3 = list(c(0.00107470281177, -93.9413043478), 
            c(0.00537351405887, -89.9413043478), 
            c(0.15394569460894, -85.9413043478), 
            c(0.24933105233145, -81.9413043478), 
            c(0.26249067859806, -77.9413043478), 
            c(0.18269947800149, -73.9413043478), 
            c(0.13218844584814, -69.9413043478), 
            c(0.01289643374128, -65.9413043478))
  d3 <- t(data.frame(d3))
  row.names(d3) <- c("12", "13", "14", "15", "16", "17", "18", "19")
  result = adjust.frequencies(frq, knownAlleles, adj=10)
  checkEquals(result$D3, d3)

  d3 = list(c(0.000424284502044, -93.9413043478), 
            c(0.002121422510221, -89.9413043478), 
            c(0.363380390341742, -85.9413043478), 
            c(0.098434004474273, -81.9413043478), 
            c(0.406233125048214, -77.9413043478), 
            c(0.072128365347528, -73.9413043478), 
            c(0.052186993751446, -69.9413043478), 
            c(0.005091414024531, -65.9413043478))
  d3 <- t(data.frame(d3))
  row.names(d3) <- c("12", "13", "14", "15", "16", "17", "18", "19")
  result = adjust.frequencies(frq, knownAlleles, adj=10, fst=0.45)
  checkEquals(result$D3, d3)
})
