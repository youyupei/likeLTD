# Case we are going to be looking at.
library(gtools)
library(microbenchmark)
library(stats)
source('R/genetics.R')
source('R/allele_report.R')
source('R/optimized.R')
source('R/scenario.R')
source('R/objectives.R')
source('R/maximize.R')
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
# defenseScenario     = do.call(defense.scenario, args)

pOpti0 = optimization.params(prosecutionScenario, nUnknowns=0, doDropin=TRUE)
pOpti1 = optimization.params(prosecutionScenario, nUnknowns=1, doDropin=TRUE)
pOpti2 = optimization.params(prosecutionScenario, nUnknowns=2, doDropin=TRUE)

# # Create objective functions for defense and prosecution.
# # prosecutionObjective = create.likelihood.vectors(prosecutionScenario, TRUE)
# defenseObjective     = create.likelihood.vectors(defenseScenario, TRUE)

arguments = list(rcont=c(0.923913043478261, 0.565217391304348,
                         0.543478260869565),
                 dropin = 0.1,
                 degradation=rep(1e-3, 4),
                 localAdjustment=1,
                 tvedebrink=-4.35,
                 dropout=c(0.175, 0.105) )
# arguments$rcont = arguments$rcont[1:(2+args$nUnknowns)]
# arguments$degradation = rep(3e-3, 3+args$nUnknowns)
# objective <- create.likelihood.vectors(defenseScenario, TRUE)

# funcs <- attr(objective, "functions") 
# a <- do.call(objective, arguments)

# bench <- microbenchmark( D3=do.call(functions$D3, arguments), 
#                          vWA=do.call(functions$vWA, arguments), 
#                          D16=do.call(functions$D16, arguments), 
#                          D2=do.call(functions$D2, arguments), 
#                          D8=do.call(functions$D8, arguments), 
#                          D21=do.call(functions$D21, arguments), 
#                          D18=do.call(functions$D18, arguments), 
#                          D19=do.call(functions$D19, arguments), 
#                          TH01=do.call(functions$TH01, arguments), 
#                          FGA=do.call(functions$FGA, arguments) )
# if (require("ggplot2")) {
#   plt <- ggplot2::qplot(y=time, data=bench, colour=expr)
#   plt <- plt + ggplot2::scale_y_log10()
#   print(plt)
# }
modify.args = function(index, newvalue, args) { 
  for(n in names(args)) {
    if(index == 1 & length(args[[n]]) == 1) { args[[n]] = newvalue; break }
    else if( index < length(args[[n]]) ) { args[[n]][index] = newvalue; break }
    else index = index - length(args[[n]])
  }
  return(args)
}

plot_all.log <- function(which=1, x=(1:99)/100.0) {
  require("ggplot2")
  
  funcme <- function(i) { 
    args$rcont = args$rcont[1:3+args$nUnknowns]
    args$degradation = args$degradation[1:3+args$nUnknowns]
    result <- do.call(prosecutionObjective, modify.args(which, x[i], args))
    result$objectives * result$penalties
  }
  y = log(sapply(1:length(x), funcme) )
  total = colSums(y)
  total = total - min(total) 
  y = rbind(y, total=total)
  data = data.frame(prob=c(y), locus=rep(rownames(y), ncol(y)),
                    rcont=rep(x, rep(nrow(y), length(x))))
  return(ggplot(data, aes(x=rcont, y=prob, group=locus, colour=locus)) + geom_line())
}

plot_all2d.log <- function(which=c(1, 2), x=(1:99)/100.0, y=(1:99)/100.0) {
  require("ggplot2")
  result = matrix(1, nrow=length(x), ncol=length(y))
  for(i in 1:length(x)) 
  {
    argi = modify.args(which[1], x[i], arguments)
    for(j in 1:length(y))
    {
      argij = modify.args(which[2], x[j], argi)
      r = do.call(objective, argij)
      result[i, j] = sum(log(r$objectives) + log(r$penalties))
    }
  }
  amap = data.frame(x=rep(x, length(y)), y=rep(y, rep(length(x), length(y))), z=c(result))
  return(ggplot(amap, aes(x, y, z=z)))
}

# arguments$degradation = log10(arguments$degradation)
# optimizeme = fixed.params(objective, arguments, zero=1e-12)

# res <- optim( optimizeme$args, optimizeme$objective, 
#               method="L-BFGS-B",
#               upper=optimizeme$upper.bounds,
#               lower=optimizeme$lower.bounds,
#               hessian=FALSE,
#               control=list(fnscale=-1, factr=1e12) )
