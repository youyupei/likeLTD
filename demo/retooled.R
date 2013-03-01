# Case we are going to be looking at.
# source('R/genetics.R')
# source('R/allele_report.R')
# source('R/optimized.R')
# source('R/scenario.R')
# source('R/objectives.R')
caseName = 'hammer'
datapath = file.path(file.path('inst', 'extdata'), caseName)

args = list(
  databaseFile = NULL,
  mixedFile    = file.path(datapath, 'hammer-CSP.csv'),
  refFile      = file.path(datapath, 'hammer-reference.csv'),
  nUnknowns    = 1,
  doDropin     = TRUE,
  ethnic       = "EA1",
  adj          = 1.0,
  fst          = 0.02,
  relatedness  = c(0, 0)
)


# Create scenarios for defense and prosecution.
# prosecutionScenario = do.call(prosecution.scenario, args)
defenseScenario     = do.call(defense.scenario, args)

# Create objective functions for defense and prosecution.
# prosecutionObjective = create.likelihood.vectors(prosecutionScenario)
# defenseObjective     = create.likelihood.vectors(defenseScenario)

arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                         1.000000000000000, 0.543478260869565),
                 dropin = 1e0, #0.108695652173913,
                 degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
                 localAdjustment=1,
                 tvedebrink=-4.35,
                 dropout=c(0.175, 0.105) )
scenarioTH01 = transform.to.locus.centric(defenseScenario)$TH01
obTH01 <- create.likelihood.per.locus(scenarioTH01)
