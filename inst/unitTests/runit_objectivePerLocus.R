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


test_empty.alleles = svTest(function() {
  
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus <- getFromNamespace("all.genotypes.per.locus",
                                                "likeLTD")
  if(! "empty.alleles" %in% ls(.GlobalEnv))
    empty.alleles <- getFromNamespace("empty.alleles", "likeLTD")
  genotypes = matrix(nrow=4, ncol=0)
  dropoutProfs = cbind(c(TRUE, FALSE, FALSE, TRUE),
                       c(FALSE, TRUE, FALSE, FALSE))
  result = empty.alleles(genotypes, dropoutProfs, 0)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 0)
  checkTrue(nrow(result) == 4)
  
  result = empty.alleles(genotypes, dropoutProfs, 1)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == 0)
  checkTrue(nrow(result) == 4)

  genotypes = all.genotypes.per.locus(4, 2)
  result = empty.alleles(genotypes, dropoutProfs, 0)
  checkTrue(is.matrix(result))
  checkTrue(ncol(result) == ncol(genotypes))
  checkTrue(nrow(result) == 4)
  checkTrue(all(result == c(FALSE, FALSE, TRUE, FALSE)))

  result = empty.alleles(genotypes, dropoutProfs, 1)
  indices = c(1, 2, 4, 5, 7, 10, 11, 12, 14, 15, 17, 20, 31, 32, 34, 35, 37,
              40, 41, 42, 44, 45, 47, 50, 61, 62, 64, 65, 67, 70, 91, 92, 94,
              95, 97, 100)
  checkTrue(ncol(result) == ncol(genotypes))
  checkTrue(nrow(result) == 4)
  checkTrue(all(!result[1, ]))
  checkTrue(all(!result[2, ]))
  checkTrue(all(!result[4, ]))
  checkTrue(all(result[3, indices]))
  result[3, indices] = FALSE
  checkTrue(all(!result))
})

test_TH01.regression.with.dropin = svTest(function() {
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
    relatedness  = c(0.5, 0)
  )
  scenario = do.call(prosecution.scenario, args)
  if(! "transform.to.locus.centric" %in% ls(.GlobalEnv))
    transform.to.locus.centric <-
      getFromNamespace("transform.to.locus.centric", "likeLTD")
  scenarioTH01 = transform.to.locus.centric(scenario)$TH01
  
  arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                           1.000000000000000, 0.543478260869565),
                   dropin = 1e0, #0.108695652173913,
                   degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                   localAdjustment=1,
                   tvedebrink=-4.35,
                   dropout=c(0.175, 0.105) )

  if(! "create.likelihood.per.locus" %in% ls(.GlobalEnv))
    create.likelihood.per.locus <-
      getFromNamespace("create.likelihood.per.locus", "likeLTD")
  scenarioTH01$nUnknowns = 2
  scenarioTH01$doDropin = TRUE
  objective.function <- create.likelihood.per.locus(scenarioTH01)

  checkEquals(do.call(objective.function, arguments), 0.0204764693571788)
# arguments$degradation = rep(2e-2, 4)
# checkEquals(do.call(objective.function, arguments), 0.0017054449886482)
# arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
#                           0.00216652022387600, 0.00131485405088635)
# checkEquals(do.call(objective.function, arguments), 0.0147496568283615)
})

test_TH01.regression.no.dropin = svTest(function() {
  # Case we are going to be looking at.
  cspPresence = matrix(c(0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0), nrow=2)
  profPresence = matrix(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=4)
  uncPresence = matrix(0, nrow=2, ncol=7)
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
  objective.function <- create.likelihood.per.locus(profPresence, cspPresence,
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
  objective.function <- create.likelihood.per.locus(profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 0, TRUE)

  checkEquals(do.call(objective.function, arguments), 9.09495530595364e-06)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 5.75530620496417e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 5.13626195539709e-05)

  objective.function <- create.likelihood.per.locus(profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 1, TRUE)
  arguments$degradation = rep(3e-3, 4)
  checkEquals(do.call(objective.function, arguments), 9.06669340994184e-06)
  arguments$degradation = rep(2e-2, 4)
  checkEquals(do.call(objective.function, arguments), 4.27972472968122e-05)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 3.77630662967064e-05)

  objective.function <- create.likelihood.per.locus(profPresence, cspPresence,
                                                    uncPresence, missingReps,
                                                    alleleDb, 2, TRUE)
  arguments$degradation = c(0.00723217060006922, 0.00569441925951047,
                            0.00216652022387600, 0.00131485405088635)
  checkEquals(do.call(objective.function, arguments), 1.3537562256385e-05)
})


test_relatedness.factors <- svTest(function() {

  alleleDb   = matrix(c(1, 0), nrow=14, ncol=2)
  row.names(alleleDb) = c("10", "11", "12", "13", "14", "14.2", "15", "16",
                          "17", "18", "19", "20", "21", "22")
  
  queriedAlleles = c("11", "21")
  alleleDb[2, 1] = 2e0
  alleleDb[13, 1] = 3e0
  
  if(! "all.genotypes.per.locus" %in% ls(.GlobalEnv))
    all.genotypes.per.locus <- getFromNamespace("all.genotypes.per.locus",
                                                "likeLTD")
  if(! "relatedness.factors" %in% ls(.GlobalEnv))
    relatedness.factors <- getFromNamespace("relatedness.factors", "likeLTD")

  genotypes = all.genotypes.per.locus(nrow(alleleDb), 2)
  # No relatedness
  result = relatedness.factors(genotypes, alleleDb, queriedAlleles, c(0, 0))
  checkEquals(result, 1)
  # First, only one relatedness 
  result = relatedness.factors(genotypes, alleleDb, queriedAlleles, c(0.9, 0))
  checkEquals(length(result), ncol(genotypes))
  na = nrow(alleleDb)
  nComb = na * (na+1) / 2
  nCont = 2
  nTrue = (2*na - 1) * nComb^(nCont-1)
  checkEquals(sum(abs(result - 1.0 * (1.0 - sum(c(0.9, 0)))) > 1e-8), nTrue)
  hasFirst = genotypes[1, ] %in% 2 | genotypes[2, ] %in% 2
  hasSecond = genotypes[1, ] %in% 13 | genotypes[2, ] %in% 13
  het = genotypes[1, ] == genotypes[2, ]

  r = result[hasFirst & (!hasSecond) & (!het)]
  check = (1-0.9) * (1e0 + 0.9 * 0.5 * 0.5 / alleleDb[2, 1])
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  r = result[hasFirst & (!hasSecond) & het]
  check = (1-0.9) * (1e0 + 0.9 * 0.5 / alleleDb[2, 1])
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  r = result[(!hasFirst) & hasSecond & (!het)]
  check = (1-0.9) * (1e0 + 0.9 * 0.5 * 0.5 / alleleDb[13, 1])
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  check = (1-0.9) * (1e0 + 0.9 * 0.5 / alleleDb[13, 1])
  r = result[(!hasFirst) & hasSecond & het]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)
  
  r = result[hasFirst & hasSecond]
  check = (1-0.9) * (1e0 + 0.9 * 0.5 * 0.5 / alleleDb[2, 1] 
                     + 0.9 * 0.5 * 0.5 / alleleDb[13, 1])
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  result = relatedness.factors(genotypes, alleleDb, queriedAlleles, c(0, 0.9))
  check = (1-0.9) * (1e0 + 0.9 * 0.5 / alleleDb[2, 1] / alleleDb[13, 1])
  r = result[hasFirst & hasSecond]
  checkTrue(all(abs(r - check) < 1e-8))
  checkTrue(length(abs(r - check) < 1e-8) > 1)

  queriedAlleles = c("2", "2")
})
