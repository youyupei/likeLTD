## Test unit 'allele_report'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests//runitallele_report.R"
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


###################################
# Finally, the unit-test themselves
###################################

"test_pack.admin.input" <- svTest(function() {
  # Tests pack.admin.input interface
  # 
  # Checks it packs the data as expected.
  # Checks it raises an exception if files do not exist and checkFile is TRUE.
  # Checks it does not raise an exception if checkFiles is FALSE.
  # Checks it does not raise an exception if checkFiles is TRUE and files exist
  databaseFile = 'data/lgc-allele-freqs-wbp.txt'
  mixedFile = 'hammer/hammer-CSP.csv'
  refFile = 'hammer/hammer-reference.csv'
  caseName = 'hammer'
  outputPath = 'hammer'

  callme <- function(checkfile) {
    # Helper function to convserve finger strength
    pack.admin.input( caseName=caseName,
                      databaseFile=databaseFile,
                      mixedFile=mixedFile,
                      refFile=refFile, 
                      outputPath=outputPath,
                      checkFiles=checkfile )
  }
 
  # Checks it packs the data as expected.
  admin <- callme(FALSE)
  checkEquals(length(admin), 5)
  checkEquals(admin$caseName, caseName)
  checkEquals(admin$databaseFile, databaseFile)
  checkEquals(admin$mixedFile, mixedFile)
  checkEquals(admin$refFile, refFile)
  checkEquals(admin$outputPath, outputPath)

  checkException(callme(checkfile=TRUE), msg="No file found exception")
  checkNoException(callme(TRUE), msg="Files not checked.")

  temporary.directory({
    dir.create('data')
    dir.create('hammer')
    file.create(file.path('data', 'lgx-allele-freqs-wbp.txt'))
    file.create(file.path('hammer', 'hammer-CSP.csv'))
    file.create(file.path('hammer', 'hammer-reference.csv'))
    checkNoException(callme(TRUE), msg="Files checked and exist.")
  })
})


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
csp.data <- function() {
  # Crime Scene Profile data used throughout the tests.

  csp = list( D3=c("14,16", "", "14,16", ""), 
               vWA=c("15,16,19", "", "15,16,17,19", ""),
               D16=c("11,13,14", "", "11,13,14", ""),
               D2=c("20,23,24,25", "", "20,24,25", ""),
               D8=c("11,12,13,15", "", "11,12,13,15", ""),
               D21=c("28,31", "", "28,29,30,31,31.2", ""),
               D18=c("", "", "13,14,16,17", ""),
               D19=c("12,13,15.2,17.2", "", "12,13,14,15.2,17.2",""),
               TH01=c("6,8,9,9.3", "", "6,8,9,9.3", ""),
               FGA=c("22", "", "22,23,25", "") )
  csp = data.frame(csp, stringsAsFactors=FALSE)
  return(csp)
}

test_read.csp.profile <- svTest(function(){
  # Tests that CSP profile is read correctly.
  
  # Path to data files.                                 
  path = Reduce(file.path, c("extdata", "hammer", "hammer-CSP.csv"))
  path = system.file(path, package="likeLTD")
  # Now read data.
  if(! "read.csp.profile" %in% ls(.GlobalEnv))
    read.csp.profile = getFromNamespace("read.csp.profile", "likeLTD")
  result = read.csp.profile(path)
  # And check it
  checkEquals(result, csp.data())
})


test_read.ref.profile <- svTest(function(){
  # Tests reading reference profiles.
  # Path to data files.                                 
  path = Reduce(file.path, c("extdata", "hammer", "hammer-reference.csv"))
  path = system.file(path, package="likeLTD")
  # Now read data.
  if(! "read.ref.profile" %in% ls(.GlobalEnv))
    read.ref.profile = getFromNamespace("read.ref.profile", "likeLTD")
  result = read.ref.profile(path)
  checkEquals(result, ref.data())
})

test_queried.vs.known <- svTest(function() {
  # Tests reading queried vs known profiles
  # Path to data files.                                 
  path = Reduce(file.path, c("extdata", "hammer", "hammer-reference.csv"))
  path = system.file(path, package="likeLTD")
  # Now read data.
  if(! "queried.vs.known" %in% ls(.GlobalEnv))
    queried.vs.known = getFromNamespace( "queried.vs.known", "likeLTD")
  result = queried.vs.known(path)
  checkEquals(result, c(TRUE, FALSE, FALSE))
})


test_internal.representation <- svTest(function() {
  # Tests that the internal representation is correctly constructed.
  if(! "internal.representation" %in% ls(.GlobalEnv))
    internal.representation = getFromNamespace("internal.representation", 
                                               "likeLTD")
  result = internal.representation(csp.data())
  checkEquals(result, internal.representation.data())
})

test_estimates.csp <- svTest(function() {
  #  Test estimate function.

  cprofs = internal.representation.data()

  check = list(a=c(85, 70, 40), b=c(100, 100, 65), c=c(92, 85, 52))
  check = data.frame(check, row.names=c("Suspect", "Victim 1", "Victim 2"))
  colnames(check) <- c("run 1", "run 2", "Total")

  if(! "estimates.csp" %in% ls(.GlobalEnv))
    estimates.csp = getFromNamespace("estimate.csp", "likeLTD")
  result = estimates.csp(ref.data(), cprofs)
  checkEquals(result, check)
})

test_summary.generator <- svTest(function() {
  # Tests summary generation unit.
  queried = c("Suspect")
  ref = ref.data()
  cprofs = internal.representation.data()
   
  # Test with known = Victim 1 only.
  if(! "summary.generator" %in% ls(.GlobalEnv))
    summary.generator = getFromNamespace( "summary.generator", "likeLTD")
  result = summary.generator(queried, c("Victim 1"), ref, cprofs)

  summary = list( D3   = c("14 16{}[]", "16 16{}[]", "{}"),
                  vWA  = c("15 19{}[]", "15 16{}[]", "{17}"),
                  D16  = c("11 14{}[]", "13 13{}[]", "{}"),
                  D2   = c("24 25{}[]", "20 20{}[]", "{23}"),
                  D8   = c("12 13{}[]", "11 15{}[]", "{}"),
                  D21  = c("28 31{}[]", "{29 30}[]", "{31.2}"),
                  D18  = c("{14 17}[]", "{17 17}[]", "{13 16}"),
                  D19  = c("15.2 17.2{}[]", "12{14}[]", "13{}"),
                  TH01 = c("9 9.3{}[]", "6 8{}[]", "{}"),
                  FGA  = c("22{23}[]", "22{25}[]", "{}") )
  summary = data.frame( summary, stringsAsFactors=FALSE,
                        row.names=c( "Suspect (Q)", "Victim 1 (K)",
                                    "Unattributable") )
  check = list( summary=summary,
                otherRep=c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
                otherUnrep=c(0, 1, 0, 1, 0, 1, 2, 0, 0, 0) )

  checkEquals(result, check)

  # Test with known = Victim 2 only.
  result = summary.generator(queried, c("Victim 2"), ref, cprofs)

  summary = list( D3   =c("14 16{}[]", "{}[15 17]", "{}"),
                  vWA  =c("15 19{}[]", "16 19{}[]", "{17}"),
                  D16  =c("11 14{}[]", "13{}[12]", "{}"),
                  D2   =c("24 25{}[]", "25{}[18]", "20{23}"),
                  D8   =c("12 13{}[]", "11 13{}[]", "15{}"),
                  D21  =c("28 31{}[]", "{29 30}[]", "{31.2}"),
                  D18  =c("{14 17}[]", "{17}[15]", "{13 16}"),
                  D19  =c("15.2 17.2{}[]", "{14 14}[]", "12 13{}"),
                  TH01 =c("9 9.3{}[]", "6{}[7]", "8{}"),
                  FGA  =c("22{23}[]", "22{}[20]", "{25}") )
  summary = data.frame( summary, stringsAsFactors=FALSE,
                        row.names=c( "Suspect (Q)", "Victim 2 (K)",
                                    "Unattributable") )
  check = list( summary=summary,
                otherRep=c(0, 0, 0, 1, 1, 0, 0, 2, 1, 0),
                otherUnrep=c(0, 1, 0, 1, 0, 1, 2, 0, 0, 1) )
  checkEquals(result, check)

  # Test with known = Victim 1&2 only.
  result = summary.generator(queried, c("Victim 1", "Victim 2"), ref, cprofs)

  summary = list( D3   = c("14 16{}[]", "16 16{}[]", "{}[15 17]", "{}"), 
                  vWA  = c("15 19{}[]", "15 16{}[]", "16 19{}[]", "{17}"), 
                  D16  = c("11 14{}[]", "13 13{}[]", "13{}[12]", "{}"), 
                  D2   = c("24 25{}[]", "20 20{}[]", "25{}[18]", "{23}"), 
                  D8   = c("12 13{}[]", "11 15{}[]", "11 13{}[]", "{}"), 
                  D21  = c("28 31{}[]", "{29 30}[]", "{29 30}[]", "{31.2}"), 
                  D18  = c("{14 17}[]", "{17 17}[]", "{17}[15]", "{13 16}"), 
                  D19  = c("15.2 17.2{}[]", "12{14}[]", "{14 14}[]", "13{}"), 
                  TH01 = c("9 9.3{}[]", "6 8{}[]", "6{}[7]", "{}"), 
                  FGA  = c("22{23}[]", "22{25}[]", "22{}[20]", "{}") )
  summary = data.frame( summary, stringsAsFactors=FALSE,
                        row.names=c( "Suspect (Q)", "Victim 1 (K)",
                                     "Victim 2 (K)", "Unattributable") )
  check = list( summary=summary,
                otherRep=c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
                otherUnrep=c(0, 1, 0, 1, 0, 1, 2, 0, 0, 0) )
  checkEquals(result, check)
})


test_unusual.alleles <- svTest(function() {
  # Tests unusual allele function
  dummyEnv <- new.env()
  data('lgc-allele-freqs-wbp', package='likeLTD', envir=dummyEnv)
  afreq      = dummyEnv[['lgc-allele-freqs-wbp']]
  
  if(! "unusual.alleles" %in% ls(.GlobalEnv))
    unusual.alleles = getFromNamespace("unusual.alleles", "likeLTD")
  result = unusual.alleles(afreq, ref.data())
  check = list( name=c("Suspect"), locus=c("D19"), allele=c("17.2"), EA1=c(0),
                EA3=c(2), EA4=c(1) )
  check = data.frame(check)
  checkEquals(result, check)
})
