library(gtools)
library(microbenchmark)
library(stats)
library(likeLTD)

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
  relatedness  = c(0, 0)/4
)

# Create scenarios for defense and prosecution.
prosecutionScenario = do.call(prosecution.scenario, args)
defenseScenario     = do.call(defense.scenario, args)

bench.any <- function(scenario, times=100L, ...) {

  objective <- create.likelihood.vectors(scenario, addAttr=TRUE, ...)
  funcs <- attr(objective, "functions") 
  arguments <- initial.arguments(scenario, ...)
  arguments$localAdjustment = 1
  microbenchmark( D3=do.call(funcs$D3, arguments), 
                  vWA=do.call(funcs$vWA, arguments), 
                  D16=do.call(funcs$D16, arguments), 
                  D2=do.call(funcs$D2, arguments), 
                  D8=do.call(funcs$D8, arguments), 
                  D21=do.call(funcs$D21, arguments), 
                  D18=do.call(funcs$D18, arguments), 
                  D19=do.call(funcs$D19, arguments), 
                  TH01=do.call(funcs$TH01, arguments), 
                  FGA=do.call(funcs$FGA, arguments), times=times )

}
bench.prosecution <- function(...) bench.any(prosecutionScenario, ...)
bench.defense <- function(...) bench.any(defenseScenario, ...)

cat(paste("Number of threads", .Call(.cpp.nbthreads, PACKAGE="likeLTD")))
cat("\n\nPROSECUTION -- unknown=0, doDropin=TRUE\n")
print(bench.prosecution(nUnknowns=0, doDropin=TRUE))
cat("\n\nPROSECUTION -- unknown=1, doDropin=TRUE\n")
print(bench.prosecution(nUnknowns=1, doDropin=TRUE))
cat("\n\nPROSECUTION -- unknown=2, doDropin=TRUE\n")
print(bench.prosecution(nUnknowns=2, doDropin=TRUE))
cat("\n\nPROSECUTION -- unknown=3, doDropin=TRUE\n")
print(bench.prosecution(nUnknowns=3, doDropin=TRUE, times=5))
