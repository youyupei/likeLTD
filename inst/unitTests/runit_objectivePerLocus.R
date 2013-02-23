## Test unit 'objectives per locus'
library(svUnit)

###############################################################
# The new two functions are to set up the unit test environment
###############################################################

.setUp <-
function () {
	## Specific actions for svUnit: prepare context
	if ("package:svUnit" %in% search()) {
		.Log <- Log() ## Make sure .Log is created
		.Log$..Unit <- "inst/unitTests/runit_objectivePerLocus.R"
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


test_TH01.regression.with.dropin = svTest(function() {
  # Case we are going to be looking at.
  cspPresence = matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0), nrow=2)
  profPresence = matrix(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4)
  queriedPresence = NULL
  uncPresence = matrix(0, nrow=2, ncol=ncol(cspPresence))
  missingReps = rep(FALSE, 2)
  alleleDb   = matrix(c( 0.00109678574626197504, 0.23251857820753871198,
                        0.17109857641686809782, 0.09871071716357775194,
                        0.15122213268869189040, 0.33767570955322767645,
                        0.00767750022383382504, -39.94130434782607608213,
                        -35.94130434782607608213, -31.94130434782607608213,
                        -27.94130434782607608213, -23.94130434782607608213,
                        -20.94130434782607608213, -19.94130434782607608213),
                      ncol=2)
  row.names(alleleDb) = c("5", "6", "7", "8", "9", "9.3", "10")

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565,
                           0.108695652173913), 
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )

  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  objective.function <- create.likelihood.per.locus(queriedPresence,
                                                    profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 2, TRUE)

  checkEquals(do.call(objective.function, arguments), 0.0204764693571788)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 0.0017054449886482)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 0.0147496568283615)
})

test_TH01.regression.no.dropin = svTest(function() {
  # Case we are going to be looking at.
  cspPresence = matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0), nrow=2)
  profPresence = matrix(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4)
  uncPresence = matrix(0, nrow=2, ncol=7)
  queriedPresence = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0), nrow=2)
  missingReps = rep(FALSE, 2)
  alleleDb   = matrix(c( 0.00109678574626197504, 0.23251857820753871198,
                        0.17109857641686809782, 0.09871071716357775194,
                        0.15122213268869189040, 0.33767570955322767645,
                        0.00767750022383382504, -39.94130434782607608213,
                        -35.94130434782607608213, -31.94130434782607608213,
                        -27.94130434782607608213, -23.94130434782607608213,
                        -20.94130434782607608213, -19.94130434782607608213),
                      ncol=2)
  row.names(alleleDb) = c("5", "6", "7", "8", "9", "9.3", "10")

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565,
                           0.108695652173913), 
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )

  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  objective.function <- create.likelihood.per.locus(queriedPresence,
                                                    profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 2, FALSE)

  checkEquals(do.call(objective.function, arguments), 0.0194512081797547)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 0.00166299033316889)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 0.0140009673609186)
})

test_D18.regression.with.dropin = svTest(function() {
  # Case we are going to be looking at.
  cspPresence = matrix(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=2)
  profPresence = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0), nrow=4)
  uncPresence = matrix(0, nrow=2, ncol=ncol(cspPresence))
  queriedPresence = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                             1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=2,)
  missingReps = rep(FALSE, 2)
  alleleDb   = matrix(c( 1.31614289551437e-02, 5.48392873130988e-03,
                         1.51356432984153e-01, 1.23936789327603e-01,
                         1.90706419554123e-01, 1.09678574626198e-03,
                         1.36001432536485e-01, 1.17356074850031e-01,
                         1.19415346047095e-01, 7.56782164920763e-02,
                         3.94842868654311e-02, 1.64517861939296e-02,
                         7.67750022383383e-03, 2.19357149252395e-03,
                         7.00586956521739e+01, 7.40586956521739e+01,
                         7.80586956521739e+01, 8.20586956521739e+01,
                         8.60586956521739e+01, 8.80586956521739e+01,
                         9.00586956521739e+01, 9.40586956521739e+01,
                         9.80586956521739e+01, 1.02058695652174e+02,
                         1.06058695652174e+02, 1.10058695652174e+02,
                         1.14058695652174e+02, 1.18058695652174e+02), ncol=2)

  row.names(alleleDb) = c("10", "11", "12", "13", "14", "14.2", "15", "16",
                          "17", "18", "19", "20", "21", "22")

  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565,
                           0.108695652173913), 
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )

  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  objective.function <- create.likelihood.per.locus(queriedPresence,
                                                    profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 0, TRUE)

  checkEquals(do.call(objective.function, arguments), 9.09495530595364e-06)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 5.75530620496417e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 5.13626195539709e-05)

  objective.function <- create.likelihood.per.locus(queriedPresence,
                                                    profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 1, TRUE)
  arguments$degradation = rep(3e-3, 4)
  checkEquals(do.call(objective.function, arguments), 9.06669340994184e-06)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 4.27972472968122e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 3.77630662967064e-05)

  objective.function <- create.likelihood.per.locus(queriedPresence,
                                                    profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 2, TRUE)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 1.3537562256385e-05)
})
