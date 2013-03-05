# Case we are going to be looking at.
source('R/genetics.R')
source('R/allele_report.R')
source('R/optimized.R')
source('R/scenario.R')
source('R/objectives.R')
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
prosecutionObjective = create.likelihood.vectors(prosecutionScenario, TRUE)
defenseObjective     = create.likelihood.vectors(defenseScenario, TRUE)

arguments = list(rcont=c(0.923913043478261, 0.565217391304348, 1.0, 
                         0.543478260869565),
                 dropin = 0.1,
                 degradation=rep(3e-3, 4),
                 localAdjustment=1,
                 tvedebrink=-4.35,
                 dropout=c(0.175, 0.105) )
objective <- create.likelihood.vectors(defenseScenario, TRUE)
funcs <- attr(objective, "functions") 
a <- do.call(objective, arguments)
# modify.args = function(index, newvalue, args) { 
#   for(n in names(args)) {
#     if(index == 1 & length(args[[n]]) == 1) { args[[n]] = newvalue; break }
#     else if( index < length(args[[n]]) ) { args[[n]][index] = newvalue; break }
#     else index = index - length(args[[n]])
#   }
#   return(args)
# }

# args = list(rcont=c(0.923913043478261, 0.565217391304348, 1.000000000000000,
#                     0.543478260869565, 0.108695652173913), 
#             dropin = 0.543478260869565, 
#             degradation=c(3e-3, 3e-3, 3e-3, 3e-3),
#             localAdjustment=1, tvedebrink=-4.35,
#             dropout=c(0.175, 0.105), beta=-4.35 )


# plot_all.log <- function(which=1, x=(1:99)/100.0) {
#   require("ggplot2")
#   
#   funcme <- function(i) { 
#     args$rcont = args$rcont[1:4]
#     args$degradation = args$degradation[1:4]
#     result <- do.call(prosecutionObjective, modify.args(which, x[i], args))
#     result$objectives * result$penalties
#   }
#   y = log(sapply(1:length(x), funcme) )
#   total = colSums(y)
#   total = total - min(total) 
#   y = rbind(y, total=total)
#   data = data.frame(prob=c(y), locus=rep(rownames(y), ncol(y)),
#                     rcont=rep(x, rep(nrow(y), length(x))))
#   return(ggplot(data, aes(x=rcont, y=prob, group=locus, colour=locus)) + geom_line())
# }
