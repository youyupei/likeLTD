estimates <- function(indiv, csp) {
  # Estimation of dropout values.
  #
  # Parameters: 
  #    indiv: Profile of an individual as a list of locus, each element
  #           containing one or two strings naming the alleles.
  #    csp: Crime-scene profile. Matrix of replicate (rows) vs loci (columns).
  #         Each element is a vector of strings naming the alleles present in
  #         the crime-scene profile. 
  meanrep = array(0, nrow(csp))
  for(rep in 1:nrow(csp))  {
    represented = c()
    for(locus in colnames(csp)) 
      if(!is.null(csp[[rep, locus]])) {
        represented = c(represented, indiv[[locus]] %in% csp[[rep, locus]])
      }
    meanrep[[rep]] <- sum(represented) / length(represented) 
  }
  meanrep <- 1e0 - meanrep / nrow(csp)
  meanrep[meanrep > 0.99] = 0.99
  meanrep[meanrep < 0.01] = 0.01
  meanrep
}

# Best(?) guess for initial arguments. 
initial.arguments <- function(hypothesis, ...) {
  # Best(?) guess for initial arguments. 
  #
  # Parameters: 
  #    hypothesis: Hypothesis for which to guess nuisance paramters. 

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  # -1 because relative to first.
  nrcont          = nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns - 1
  rcont           = rep(1, nrcont)
  degradation     = rep(3e-3, nrcont + 1)
  localAdjustment = rep(1, ncol(hypothesis$dropoutProfs))
  dropout         = rep(0.5, nrow(hypothesis$cspProfile))

  list(rcont           = rcont,
       degradation     = degradation,
       localAdjustment = localAdjustment,
       tvedebrink      = -4.35,
       dropin          = 1e-2,
       dropout         = dropout)
}


upper.bounds = function(arguments, zero=1e-8) { 
  # Upper bounds of optimization function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimization function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  degradation = rep(Inf, length(arguments$degradation))
  rcont       = rep(Inf, length(arguments$rcont))
  localAdjustment = rep(Inf, length(arguments$localAdjustment))
  dropout     = rep(1-zero, length(arguments$dropout))

  list(rcont           = rcont, 
       dropin          = Inf,
       degradation     = degradation,
       localAdjustment = localAdjustment,
       tvedebrink      = 1-zero, 
       dropout         = dropout)[names(arguments)]
}
lower.bounds = function(arguments, zero=1e-8, logDegradation=FALSE) { 
  # Lower bounds of optimization function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimization function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  #   logDegradation: Wether degradation parameters are entered as exponents of
  #                   10.
  degradation = if(logDegradation) { -Inf } else { 0 }
  degradation = rep(degradation, length(arguments$degradation))
  rcont       = rep(zero, length(arguments$rcont))
  localAdjustment = rep(zero, length(arguments$localAdjustment))
  dropout     = rep(zero, length(arguments$dropout))

  list(rcont           = rcont, 
       dropin          = zero,
       degradation     = degradation,
       localAdjustment = localAdjustment,
       tvedebrink      = -Inf,
       dropout         = dropout)[names(arguments)]
}


optimization.params <- function(hypothesis, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=1e-8, ...) {
  # Creates the optimization parameters for optim.
  #
  # optim is the optimization function from R's stat package.
  # 
  # Parameters:
  #    hypothesis: Hypothesis for which to create the objective function.
  #    verbose: Wether to print each time the likelihood is computed.
  #    fixed: List of arguments to keep fixed.
  #    logObjective: Whether to optimize the log of the likelihood.
  #    logDegradation: Whether to input the degradation as 10^d or not.
  #    arguments: Starting arguments. If NULL, then uses initial.arguments to
  #               compute them.
  #    zero: An epsilonic number used to indicate lower and upper bounds which
  #          should be excluded.

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  objective = create.likelihood.vectors(hypothesis, ...)

  if(is.null(arguments)) arguments = initial.arguments(hypothesis)
  if(logDegradation) arguments$degradation = log10(arguments$degradation)
  
  template = arguments

  if(!is.null(fixed)) {
    fixedstuff = template[fixed]
    template = template[-which(names(arguments) %in% fixed)]
  } else  fixedstuff = NULL

  result.function <- function(x) {
 
    # If a flat vector, reconstruct list. 
    if(typeof(x) == "double") x = relist(x, template)
    # Make sure it contains fixed terms.
    if(length(fixedstuff)) x = append(x, fixedstuff)
    # Converts degradation from log10 format.
    if(logDegradation) x$degradation = 10^x$degradation
    # Compute objective function.
    result <- do.call(objective, x)
    # Compute as log if requested, otherwise as product.
    if(logObjective)
      result <- sum(log(result$objectives) + log(result$penalties))
    else result <- log(prod(result$objectives) * prod(result$penalties))
    # Print out if requested.
    if(verbose) {
      # print(unlist(append(x, list(result=result))))
      print(result)
    }
    # return result
    result
  }
  
  lower = lower.bounds(arguments, zero, logDegradation)
  upper = upper.bounds(arguments, zero)
  lower = lower[names(template)] 
  upper = upper[names(template)] 

  list(par     = unlist(template), 
       fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       control = list(fnscale=-1, factr=1e12), 
       method  = "L-BFGS-B",
       hessian = FALSE )
}
