library('stats')
library('likeLTD')


# Case we are going to be looking at.
caseName = 'hammer'
datapath = file.path(system.file("extdata", package="likeLTD"), caseName)

args = list(
  databaseFile = NULL,
  mixedFile    = file.path(datapath, 'hammer-CSP.csv'),
  refFile      = file.path(datapath, 'hammer-reference.csv'),
  nUnknowns    = 0,
  doDropin     = TRUE,
  ethnic       = "EA1",
  adj          = 1.0,
  fst          = 0.02,
  relatedness  = c(0, 0)/4
)

# Create hypothesis for defense and prosecution.
prosecutionHyp = do.call(prosecution.hypothesis, args)
defenseHyp     = do.call(defense.hypothesis, args)

# Create optimization parameters, with some modification to hypothesis
# arguments. nUnknowns is modified here and now. It does not need to be, but it
# can be. The somethingParams values contain the objective functions and
# optimization parameters.
prosecutionParams <- optimization.params(prosecutionHyp, verbose=TRUE,
                                         nUnknowns=1) 
defenseParams <- optimization.params(defenseHyp, verbose=FALSE, nUnknowns=2)


# Now perform actual optimization.
prosecutionResult <- do.call(optim, prosecutionParams)
defenseResult     <- do.call(optim, defenseParams)
