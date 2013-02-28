# Case we are going to be looking at.
caseName = 'hammer'

args = list(
  caseName     = caseName,
  datapath     = file.path(file.path('inst', 'extdata'), caseName),
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


# Create scenarios for defense and prosecution.
prosecutionScenario = do.call(prosecution.scenario, args)
defenseScenario     = do.call(defense.scenario, args)

# Create objective functions for defense and prosecution.
prosecutionObjective = create.likelihood.vectors(prosecutionScenario)
defenseObjective     = create.likelihood.vectors(defenseScenario)
