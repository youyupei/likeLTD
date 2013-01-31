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


##############################################
# Then two functions to help design unit tests
##############################################

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

###################################
# Finally, the unit-test themselves
###################################
test_known.alleles <- svTest(function()  {
  # Tests known alleles
 
  if(! "known.alleles" %in% ls(.GlobalEnv))
    known.alleles = getFromNamespace("known.alleles", "likeLTD")
  result <- known.alleles("Suspect", "Victim 1", ref.data()) 
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
  checkEquals(result, check)

  result <- known.alleles("Suspect", c("Victim 1", "Victim 2"), ref.data())
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
  checkEquals(result, check)
})

test_check.dropouts <- svTest(function() {
  # Checks that we can figure out dropouts.
  cprofs = list( D3=  list( list(csp=c("14", "16"), unc=""),
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

  if(! "has.dropouts" %in% ls(.GlobalEnv))
    has.dropouts = getFromNamespace("has.dropouts", "likeLTD")
  checkTrue(has.dropouts("Suspect", ref.data(), cprofs))

  # D18 rep 1 is empty, so this checks that limit case as well.
  cprofs$FGA[[1]]$csp <- c("22", "23")
  checkTrue(has.dropouts("Suspect", ref.data(), cprofs))
})
