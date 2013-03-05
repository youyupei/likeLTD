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
  relatedness  = c(1, 1)/4
)


# Create scenarios for defense and prosecution.
prosecutionScenario = do.call(prosecution.scenario, args)
defenseScenario     = do.call(defense.scenario, args)

# Create objective functions for defense and prosecution.
# prosecutionObjective = create.likelihood.vectors(prosecutionScenario)
# defenseObjective     = create.likelihood.vectors(defenseScenario)

arguments = list(rcont=c(1.0, 0.923913043478261, 0.565217391304348,
                         0.543478260869565),
                 dropin = 0.108695652173913,
                 degradation=rep(3e-3, 4),
                 localAdjustment=1,
                 tvedebrink=-4.35,
                 dropout=c(0.15, 0.01) )
objective <- create.likelihood.vectors(prosecutionScenario, TRUE)
funcs <- attr(objective, "functions") 
a <- do.call(objective, arguments)

# args$rcont = c(1, 0)/4
# defenseScenario     = do.call(defense.scenario, args)
arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0, 
                         0.543478260869565),
                 dropin = 0.108695652173913,
                 degradation=rep(3e-3, 4),
                 localAdjustment=1,
                 tvedebrink=-4.35,
                 dropout=c(0.175, 0.105) )
objective <- create.likelihood.vectors(defenseScenario, TRUE)
funcs <- attr(objective, "functions") 
a <- do.call(objective, arguments)
