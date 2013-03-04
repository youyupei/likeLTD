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
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    mixedFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 0,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0)
  )
  if(! "defense.scenario" %in% ls(.GlobalEnv))
    defense.scenario <- getFromNamespace("defense.scenario", "likeLTD")
  scenario = do.call(defense.scenario, args)
  scenario$nUnknowns = 0

  likelihood <- create.likelihood.vectors(scenario)
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

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348),
                   dropin = 1.0,
                   degradation=rep(3e-3, 2),
                   localAdjustment=1,
                   tvedebrink=-4.35, beta=-4.35,
                   dropout=c(0.175, 0.105) )
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, check$objectives)
})
