library('ggplot2')
library('stats')

# Case we are going to be looking at.
caseName = 'hammer'
datapath = file.path(system.file("extdata", package="likeLTD"), caseName)

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

# Create scenarios for defense
defenseScenario = do.call(defense.scenario, args)

# Import private functions.
upper.bounds <- getFromNamespace("upper.bounds", "likeLTD")
lower.bounds <- getFromNamespace("lower.bounds", "likeLTD")

plotme2d <- function(scenario, which=c(1, 2), large=100, N=20, arguments=NULL,
                     x=NULL, y=NULL) {

  # Create objective function
  objective <- create.likelihood.log(scenario, verbose=FALSE)
  # This gives us the shape of the arguments.
  if(is.null(arguments)) arguments <- initial.arguments(scenario)
  # Figure out the bounds. 
  upper <- upper.bounds(arguments)
  lower <- lower.bounds(arguments)

  # Now figure out extents.
  xRange = c(unlist(lower)[which[[1]]], unlist(upper)[which[[1]]]) 
  yRange = c(unlist(lower)[which[[2]]], unlist(upper)[which[[2]]]) 

  # Replace infinities with large number.
  if(length(large) == 1) large = rep(large, 2)
  if(abs(xRange[[1]]) == Inf) xRange[[1]] = sign(xRange[[1]]) * large[[1]]
  if(abs(xRange[[2]]) == Inf) xRange[[2]] = sign(xRange[[2]]) * large[[1]]
  if(abs(yRange[[1]]) == Inf) yRange[[1]] = sign(yRange[[1]]) * large[[2]]
  if(abs(yRange[[2]]) == Inf) yRange[[2]] = sign(yRange[[2]]) * large[[2]]
  
  # Now create map
  if(length(N) == 1) N = rep(N, 2)
  if(is.null(x))
    x = ((1:N[[1]]) - 1) / (N[[1]]-1) * (xRange[[2]] - xRange[[1]]) + xRange[[1]]
  if(is.null(y)) 
    y = ((1:N[[2]]) - 1) / (N[[2]]-1) * (yRange[[2]] - yRange[[1]]) + yRange[[1]]
  

  map = matrix(0, nrow=length(x), ncol=length(y))
  for(i in 1:length(x)) for(j in 1:length(y)) {
    newargs = unlist(arguments)
    newargs[[which[[1]]]] = x[[i]]
    newargs[[which[[2]]]] = y[[j]]
    map[i, j] = do.call(objective, relist(newargs, arguments))
  }

  amap = data.frame(x=rep(x, length(y)), y=rep(y, rep(length(x), length(y))),
                    z=c(map))
  title = sprintf("nUnknowns=%d, dropin=%d", scenario$nUnknowns,
                  scenario$doDropin)
  ggplot(amap, aes(x, y, z=z))                      +
     xlab(names(unlist(arguments))[[which[[1]]]])   +
     ylab(names(unlist(arguments))[[which[[2]]]])   +
     geom_tile(aes(fill=z))                         +
     stat_contour()                                 +
     labs(fill="Log-likelihood", title=title)
}

# plot_all2d.log(which=c(1,3), nUnknowns=1, x=(1:20/10.0), y=(1:20/10.0)) + 
#          geom_tile(aes(fill=z))                                         + 
#          stat_contour()
