## Test unit 'objectives'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_objective.R"
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


test_regression <- svTest(function() {
  # Case we are going to be looking at.
  caseName = 'hammer'
  packagepath = system.file(file.path('extdata', caseName), package="likeLTD")
  # Construct input: frequency file.
  databaseFile = NULL #file.path(datapath, 'lgc-allele-freqs-wbp.txt')
  # Construct input: crime scene profile
  mixedFile = file.path(packagepath, 'hammer-CSP.csv')
  # Construct input: reference profiles
  refFile = file.path(packagepath, 'hammer-reference.csv')
  # Construct input: output path in the R temp directory for now

  # Construct list of all administrative input
  admin <- pack.admin.input( caseName=caseName,
                             databaseFile=databaseFile,
                             mixedFile=mixedFile,
                             refFile=refFile )

  if(! "create.likelihood.vectors" %in% ls(.GlobalEnv))
    create.likelihood.vectors <-
      getFromNamespace("create.likelihood.vectors", "likeLTD")
  likelihood <- create.likelihood.vectors(admin, nUnknowns=0, doDropin=TRUE)
  objectives = c(3.10250372325746e-04, 1.17224578453062e-02,
                 4.76863464507366e-05, 4.29822197384531e-06,
                 1.76350988087800e-03, 1.43444739873622e-08,
                 9.09495530595362e-06, 4.74945076009649e-11,
                 3.01744911691531e-04, 6.89041414610852e-03)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  penalties = c(0.00303192703176332, 0.00303192703176332, 0.00303192703176332,
                0.00303192703176332, 0.00303192703176332, 0.00303192703176332, 
                0.00303192703176332, 0.00303192703176332, 0.00303192703176332,
                0.00303192703176332)

  check = list(objectives=objectives, penalties=penalties)

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565,
                           0.108695652173913), 
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35, beta=-4.35,
                   dropout=c(0.175, 0.105) )
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, check$objectives)
})
