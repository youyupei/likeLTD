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
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  # -1 because relative to first.
  nrcont          = max(nrow(hypothesis$dropoutProfs)
                        + hypothesis$nUnknowns - 1, 0)
  localAdjustment = rep(1, ncol(hypothesis$dropoutProfs))
  dropout         = runif(nrow(hypothesis$cspProfile),min=0.3,max=0.7)
  degradation     = rep( 3e-3, 
                         nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns )
  rcont           = runif(nrcont, min=0.5, max=1.5)
  dropin          = NULL
  if(hypothesis$doDropin) dropin = 1e-2

  list(localAdjustment = localAdjustment,
       tvedebrink      = -4.35,
       dropout         = dropout, 
       degradation     = degradation,
       rcont           = rcont,
       dropin          = dropin)
}

relistArguments <- function( parameters, hypothesis, fixed=NULL,
                             logDegradation=TRUE, arguments=NULL ) {
  # Remakes arguments from a flat vector into a list.
  if(is.null(arguments)) arguments = initial.arguments(hypothesis) 
  if(!is.null(fixed)) 
       template = arguments[-which(names(arguments) %in% fixed)]
  else template = arguments
  notempty = Filter(function(n) length(n) > 0, template)
  result <- relist(parameters, notempty)
  if(logDegradation && "degradation" %in% names(result))
    result[["degradation"]] = 10^result[["degradation"]]
  if(!is.null(fixed)) result <- append(result, arguments[fixed])
  result
}



upper.bounds = function(arguments, zero=1e-4) { 
  # Upper bounds of optimization function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimization function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  localAdjustment = rep(Inf, length(arguments$localAdjustment))
  dropout     = rep(1-zero, length(arguments$dropout))
  degradation = rep(Inf, length(arguments$degradation))
  rcont       = rep(Inf, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = Inf

  list(localAdjustment = localAdjustment,
       tvedebrink      = 0-zero, 
       dropout         = dropout,
       degradation     = degradation,
       rcont           = rcont,
       dropin          = dropin)[names(arguments)]
}
lower.bounds = function(arguments, zero=1e-4, logDegradation=FALSE) { 
  # Lower bounds of optimization function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimization function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  #   logDegradation: Wether degradation parameters are entered as exponents of
  #                   10.
  localAdjustment = rep(zero, length(arguments$localAdjustment))
  degradation = if(logDegradation) { -Inf } else { 0 }
  degradation = rep(degradation, length(arguments$degradation))
  dropout     = rep(zero, length(arguments$dropout))
  rcont       = rep(zero, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = zero

  list(localAdjustment = localAdjustment,
       tvedebrink      = -Inf,
       degradation     = degradation,
       dropout         = dropout,
       rcont           = rcont, 
       dropin          = dropin)[names(arguments)]
}


optimization.params <- function(hypothesis, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=1e-4, throwError=TRUE,
                                withPenalties=TRUE, ...) {
  # Creates the optimization parameters for optim.
  #
  # optim is the optimization function from R's stat package.
  # 
  # Parameters:
  #    hypothesis: Hypothesis for which to create the objective function.
  #    verbose: Wether to print each time the likelihood is computed.
  #    fixed: List of arguments to keep fixed.
  #    logObjective: Whether to optimize the log10 of the likelihood.
  #    logDegradation: Whether to input the degradation as 10^d or not.
  #    arguments: Starting arguments. If NULL, then uses initial.arguments to
  #               compute them.
  #    zero: An epsilonic number used to indicate lower and upper bounds which
  #          should be excluded.
  #    throwError: Throw an error if result is infinite

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  sanity.check(hypothesis) # makes sure hypothesis has right type.
  objective = create.likelihood.vectors(hypothesis, ...)


  args = arguments
  if(is.null(args)) args = initial.arguments(hypothesis)
  if(logDegradation && "degradation" %in% names(args)) 
    args$degradation = log10(args$degradation)
  
  # Make sure we don't include empty stuff (like rcont sometimes)
  template = Filter(function(n) length(n) > 0, args)

  if(!is.null(fixed)) {
    fixedstuff = args[fixed]
    template = args[-which(names(args) %in% fixed)]
  } else  fixedstuff = NULL

  result.function <- function(x) {
 
    # If a flat vector, reconstruct list. 
    if(typeof(x) == "double")
      x = relistArguments(x, hypothesis, fixed=fixed, arguments=arguments)
    # Otherwise, checks some options.
    else { 
      # Make sure it contains fixed terms.
      if(length(fixedstuff)) x = append(x, fixedstuff)
      # Converts degradation from log10 format.
      if(logDegradation && "degradation" %in% names(x))
        x$degradation = 10^x$degradation
    }
    # Compute objective function.
    result <- do.call(objective, x)
    # Compute as log if requested, otherwise as product.
    if(withPenalties) {
      if(logObjective)
        result <- sum(log10(result$objectives) + log10(result$penalties))
      else result <- prod(result$objectives) * prod(result$penalties)
    } else if(logObjective) {
      result <- sum(log10(result$objectives))
    } else result <- prod(result$objectives) 
    # Print out if requested.
    if(verbose) {
      # print(unlist(append(x, list(result=result))))
      print(result)
    }
    # If result is infinite, do throw
    if(throwError == TRUE && is.infinite(result)) {
      cat("Objective function is over/underflow: ", result, "\n")
      print(x)
      stop("Objective function over/underflow")
    }
    # return result
    result
  }
  
  lower = lower.bounds(args, zero, logDegradation)
  upper = upper.bounds(args, zero)
  lower = lower[names(template)] 
  upper = upper[names(template)] 


  list(par     = unlist(template), 
       fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       control = list(fnscale=-1, factr=1e7, maxit=500), 
       method  = "L-BFGS-B",
       hessian = FALSE )
}

multiRun <- function(hypothesis, nrun, ...) 
  {
	# Runs the optimisation for a given hypothesis nrun times
	# before returning the results that gave the maximum likelihood,
	# as well as the standard deviation of the likelihoods, and the
	# number of runs that completed
	results <- list(); L=NULL; startParams <- list()
	for(i in 1:nrun) 
		{
		# If the previous run did not converge, do it again
		if(i>1)
			{
			condition = results[[i-1]]$convergence!=0
			if(condition) i = i - 1
			}
		# Set parameters (including randomisation)
		params = optimization.params(hypothesis, ...)
		startParams[[i]] <- params$par
		results[[i]] <- try(do.call(optim, params), TRUE)
		# If optim returns an error record likelihood as NA
		if((class(results[[i]])!="try.error")&(length(results[[i]])!=1)) 
			{
			L <- c(L, results[[i]]$value)
			} else {
			L <- c(L, NA)
			}
		}
	# Index of which runs returned finite L
	finiteIndex <- which(is.finite(L))
	# Index of which run returned highest L
    	maxIndex <- finiteIndex[which.max(L[finiteIndex])]
    	# Produce results to be returned
    	out = list()
    	out$startParams = startParams[[maxIndex]]
	out$likeStanDev = sd(L[finiteIndex],na.rm=TRUE)
	out$numComplete = length(finiteIndex)
	out = append(out, results[[maxIndex]])
	return(out)
	}
