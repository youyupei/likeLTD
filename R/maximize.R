estimates <- function(indiv, csp) {
  # Estimation of dropout values.
  meanrep = array(0, nrow(csp))
  for(rep in 1:nrow(csp))  {
    nloci = 0e0
    for(locus in colnames(csp)) 
      if(!is.null(csp[[rep, locus]])) {
        represented = intersect(csp[[rep, locus]], indiv[[locus]])
        meanrep[[rep]] <- meanrep[[rep]] + 
          length(represented) / length(indiv[[locus]])
        nloci <- nloci + 1
      }
    meanrep[[rep]] <- meanrep[[rep]] / nloci
  }
  meanrep <- 1e0 - meanrep / nrow(csp)
  meanrep[meanrep > 0.99] = 0.99
  meanrep[meanrep < 0.01] = 0.01
}

initial.arguments <- function(scenario) {
  # Best guess for initial arguments. 
  rcont = apply(scenario$dropoutProfiles, 1,
                function(n) estimates(n, scenario$cspProfile))
  rcont <- c(rcont, rep(mean(rcont), scenario$nUnknowns))
  rcont <- rcont / rcont[[1]]
  list(dropout = estimates(scenario$queriedProfile,
                           scenario$cspProfile),
       rcont = rcont[-1],
       degradation = rep(3e-3, nrow(scenario$dropoutProfiles)),
       localAdjustment = rep(1, ncol(scenario$dropoutProfiles)),
       tvedebrink = -4.35,
       dropin = 0e0
  )
}


optimization.params <- function(scenario, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, deglog="log",
                                arguments=NULL, zero=1e-8, ...) {
  # Creates the optimization parameters for optim.
  #
  # optim is the optimization function from R's stat package.
  # 
  # Parameters:
  #    scenario: Scenario for which to create the objective function.
  scenario = add.args.to.scenario(scenario, ...)
  objective = create.likelihood.vectors(scenario, ...)

  if(is.null(arguments)) arguments = initial.arguments(scenario)
  if(deglog) arguments$degradation = log10(arguments$degradation)
  
  template = arguments

  if(!is.null(fixed)) {
    fixedstuff = template[fixed]
    notfixed = template[-which(names(arguments) %in% fixed)]
    template = notfixed
  }

  result.function <- function(x) {
 
    # If a flat vector, reconstruct list. 
    if(typeof(x) == "double") x = relist(x, template)
    # Make sure it contains fixed terms.
    if(length(fixedstuff)) x = append(x, fixedstuff)
    # Converts degradation from log10 format.
    if(deglog) x$degradation = 10^x$degradation
    # Compute objective function.
    result <- do.call(objective, x)
    # Compute as log if requested, otherwise as product.
    if(logObjective)
      result <- sum(log(result$objectives) + log(result$penalties))
    else result <- log(prod(result$objectives) * prod(result$penalties))
    # Print out if requested.
    if(verbose) print(unlist(append(x, list=(result=result))))
    # return resul
    result
  }
  return(list(par=unlist(notfixed), fn=result.function, 
              lower.bounds=lower.bounds(arguments, zero=zero, fixed=fixed),
              upper.bounds=upper.bounds(arguments, zero=zero, fixed=fixed)),
              control=list(fnscale=-1, factr=1e12), method="L-BFGS-B",
              hessian=FALSE )
}

upper.bounds = function(arguments, zero=1e-8, fixed=NULL) { 
  result = list(rcont=rep(Inf, length(arguments$rcont)),
                dropin=Inf,
                degradation=rep(Inf, length(arguments$degradation)),
                localAdjustment=Inf,
                tvedebrink=1-zero,
                dropout=rep(Inf, length(arguments$dropout)))
  if(!is.null(fixed)) result = result[-which(names(arguments) %in% fixed)]
  unlist(result)
}
lower.bounds = function(arguments, zero=1e-8, fixed=NULL) { 
  result = list(rcont=rep(zero, length(arguments$rcont)),
                dropin=zero,
                degradation=rep(-Inf, length(arguments$degradation)),
                localAdjustment=zero,
                tvedebrink=-Inf,
                dropout=rep(zero, length(arguments$dropout)))
  if(!is.null(fixed)) result = result[-which(names(arguments) %in% fixed)]
  unlist(result)
}

