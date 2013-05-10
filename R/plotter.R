# Plots two dimensional graph of likelihood
plotLikelihood.2d <- function(hypothesis, which=c(1, 2), large=100, N=20,
                              arguments=NULL, x=NULL, y=NULL,
                              logObjective=TRUE, logDegradation=TRUE,
                              contours=list()) {

  if(!require(ggplot2)) stop("Plotting reqires ggplot2")
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  # Create objective function
  creator <- create.likelihood
  if(logObjective) creator <- create.likelihood.log

  if(logDegradation) {
    objective0 <- creator(hypothesis, verbose=FALSE)
    objective <- function(...) {
      args = list(...)
      args$degradation <- 10^args$degradation
      do.call(objective0, args)
    }
  } else objective <- creator(hypothesis, verbose=FALSE)
  
  # This gives us the shape of the arguments.
  if(is.null(arguments)) {
    arguments <- initial.arguments(hypothesis)
    if(logDegradation) arguments$degradation = log10(arguments$degradation)
  }
  # Figure out the bounds. 
  upper <- upper.bounds(arguments)
  lower <- lower.bounds(arguments, logDegradation=TRUE)

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
  z = amap$z # Avoids buggy warning when checking package.
  title = sprintf("nUnknowns=%d, dropin=%d", hypothesis$nUnknowns,
                  hypothesis$doDropin)
  fill = "Likelihood"
  if(logObjective) fill = "Log-likelihood"
  ggplot(amap, aes(x=x, y=y, z=z))                  +
     xlab(names(unlist(arguments))[[which[[1]]]])   +
     ylab(names(unlist(arguments))[[which[[2]]]])   +
     geom_tile(aes(fill=z))                         +
     labs(fill=fill)                                + 
     do.call(stat_contour, contours)
}
