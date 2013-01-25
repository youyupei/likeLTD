## Test unit 'allele_report'

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
if (!"package:svUnit" %in% search()) svTest <- function(expr) { return(expr); }

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

"test_pack.admin.input" <- svTest(function() {
  # Tests pack.admin.input interface
  # 
  # Checks it packs the data as expected.
  # Checks it raises an exception if files do not exist and checkFile is TRUE.
  # Checks it does not raise an exception if checkFiles is FALSE.
  # Checks it does not raise an exception if checkFiles is TRUE and files exist
  frequencyFile = 'data/lgc-allele-freqs-wbp.txt'
  mixedFile = 'hammer/hammer-CSP.csv'
  refFile = 'hammer/hammer-reference.csv'
  caseName = 'hammer'
  outputPath = 'hammer'

  callme <- function(checkfile) {
    # Helper function to convserve finger strength
    pack.admin.input( caseName=caseName,
                      frequencyFile=frequencyFile,
                      mixedFile=mixedFile,
                      refFile=refFile, 
                      outputPath=outputPath,
                      checkFiles=checkfile )
  }
 
  # Checks it packs the data as expected.
  admin <- callme(FALSE)
  checkEquals(length(admin), 5)
  checkEquals(admin$caseName, caseName)
  checkEquals(admin$frequencyFile, frequencyFile)
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
