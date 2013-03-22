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


test_regression1 <- svTest(function() {
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
  if(! "defense.hypothesis" %in% ls(.GlobalEnv))
    defense.hypothesis <- getFromNamespace("defense.hypothesis", "likeLTD")
  hypothesis = do.call(defense.hypothesis, args)
  hypothesis$nUnknowns = 0

  likelihood <- create.likelihood.vectors(hypothesis)
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

test_regression.zerounknown <- svTest(function() {
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
  defenseHypothesis = do.call(defense.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(1.0, 0.923913043478261, 0.565217391304348),
                   dropin = 0.543478260869565,
                   degradation=c(3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(2.64877312615097e-04, 8.30010229904578e-02,
                 6.57219727208952e-02, 3.55156330421480e-03,
                 5.58990000167620e-01, 2.58579653708631e-05,
                 6.88667211400942e-07, 1.63060515953163e-05,
                 2.70886382744855e-02, 3.64303077108431e-03)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 5.82387768745086e-25)

  likelihood <- create.likelihood.vectors(defenseHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0),
                   dropin = 0.543478260869565,
                   degradation=c(3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(2.01470581026018e-03, 1.13394042752114e-02, 4.69287546906059e-03,
                 2.47293709587760e-04, 5.90787709164540e-02, 9.97729081224379e-07,
                 6.25814358291258e-06, 4.76982275758800e-08, 1.46838789834590e-02,
                 1.43741838557565e-03)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 6.97617046679607e-32)
})

test_regression.oneunknown <- svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    mixedFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 1,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(0.0, 0)
  )
  defenseHypothesis = do.call(defense.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(1.0, 0.923913043478261, 0.565217391304348,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(5.24609987450467e-05, 1.28369768074416e-01, 3.19566044368621e-02, 
                 9.37696749226556e-04, 3.54990168257584e-01, 3.61714632221897e-05, 
                 2.78024377866913e-06, 2.04906731738543e-04, 2.15249772203121e-02, 
                 4.27467268573391e-04)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 1.9762409429182e-25)

  likelihood <- create.likelihood.vectors(defenseHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(5.40928563410066e-04, 1.24028883413425e-02, 4.32400562747164e-03,
                 3.71369198539134e-04, 7.20410130310739e-02, 1.08896493758788e-06,
                 2.81245481065190e-06, 9.38982593842111e-07, 1.94351600773495e-02,
                 5.75805705402794e-04)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 3.63417836957934e-31)
})

test_regression.relatedness <- svTest(function() {
  # Case we are going to be looking at.
  datapath     = system.file(file.path('extdata', 'hammer'), package="likeLTD")
  args = list(
    databaseFile = NULL,
    mixedFile    = file.path(datapath, 'hammer-CSP.csv'),
    refFile      = file.path(datapath, 'hammer-reference.csv'),
    nUnknowns    = 1,
    doDropin     = TRUE,
    ethnic       = "EA1",
    adj          = 1.0,
    fst          = 0.02,
    relatedness  = c(1, 1)/4
  )
  defenseHypothesis = do.call(defense.hypothesis, args)
  prosecutionHypothesis = do.call(prosecution.hypothesis, args)

  likelihood <- create.likelihood.vectors(prosecutionHypothesis)
  arguments = list(rcont=c(1.0, 0.923913043478261, 0.565217391304348,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.15, 0.01) )
  objectives = c(5.24609987450467e-05, 1.28369768074416e-01, 3.19566044368621e-02, 
                 9.37696749226556e-04, 3.54990168257584e-01, 3.61714632221897e-05, 
                 2.78024377866913e-06, 2.04906731738543e-04, 2.15249772203121e-02, 
                 4.27467268573391e-04)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 1.9762409429182e-25)

  likelihood <- create.likelihood.vectors(defenseHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(1.80216091977255e-03, 3.60032286625716e-02, 3.91628936030574e-02,
                 1.35921436545295e-03, 1.70921706383154e-01, 6.80165107428541e-06,
                 2.21989141657041e-06, 4.33073720864727e-05, 5.01713248581239e-02,
                 2.04269633599325e-03)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 5.75589320779353e-25)

  defenseHypothesis$relatedness = c(0.5, 0.75)/4
  likelihood <- create.likelihood.vectors(defenseHypothesis)
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0,
                           0.543478260869565),
                   dropin = 0.108695652173913,
                   degradation=rep(3e-3, 4),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )
  objectives = c(1.42928505786127e-03, 2.89308741853859e-02, 2.91781768719226e-02,
                 1.07691206144212e-03, 1.40960034434433e-01, 5.21398235786849e-06,
                 2.38633107944920e-06, 3.20697394315575e-05, 4.11927291159638e-02,
                 1.61542861999956e-03)
  names(objectives) = c("D3", "vWA", "D16", "D2", "D8", "D21", "D18", "D19",
                        "TH01", "FGA")
  result <- do.call(likelihood, arguments)
  checkEquals(result$objectives, objectives)
  checkEquals(prod(result$objectives * result$penalties), 7.07569814196917e-26)
})
