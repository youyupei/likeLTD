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
  locusAdjustment = rep(1, ncol(hypothesis$dropoutProfs))
  dropout         = runif(nrow(hypothesis$cspProfile),min=0.3,max=0.7)
  degradation     = rep( 3e-3, 
                         nrow(hypothesis$dropoutProfs) + hypothesis$nUnknowns )
  rcont           = runif(nrcont, min=0.5, max=1.5)
  dropin          = NULL
  if(hypothesis$doDropin) dropin = 1e-2

  list(locusAdjustment = locusAdjustment,
       power           = -4.35,
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



upper.bounds = function(arguments, zero=1e-6) { 
  # Upper bounds of optimisation function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimisation function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  locusAdjustment = rep(1.5, length(arguments$locusAdjustment))
  dropout     = rep(1-zero, length(arguments$dropout))
  degradation = rep(0-zero, length(arguments$degradation))
  rcont       = rep(100, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = 10 - zero

  list(locusAdjustment = locusAdjustment,
       power           = -2, 
       dropout         = dropout,
       degradation     = degradation,
       rcont           = rcont,
       dropin          = dropin)[names(arguments)]
}
lower.bounds = function(arguments, zero=1e-6, logDegradation=FALSE) { 
  # Lower bounds of optimisation function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimisation function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  #   logDegradation: Wether degradation parameters are entered as exponents of
  #                   10.
  locusAdjustment = rep(0.5, length(arguments$locusAdjustment))
  degradation = if(logDegradation) { -20 } else { 0 }
  degradation = rep(degradation, length(arguments$degradation))
  dropout     = rep(0, length(arguments$dropout))
  rcont       = rep(zero, length(arguments$rcont))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = zero

  list(locusAdjustment = locusAdjustment,
       power           = -6,
       degradation     = degradation,
       dropout         = dropout,
       rcont           = rcont, 
       dropin          = dropin)[names(arguments)]
}


optimisation.params <- function(hypothesis, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=0, throwError=FALSE,
                                withPenalties=TRUE, objective=NULL, iterMax=NULL,...) {
  # Creates the optimisation parameters for optim.
  #
  # optim is the optimisation function from R's stat package.
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
  # If the objective function has not been handed to optimizatio.params,
  # make the objective function
  if(is.null(objective)) objective = create.likelihood.vectors(hypothesis, ...)
  
  # Load DEoptim
  if(!require(DEoptim)) install.packages("DEoptim")

  # Get maximum allele fraction
  maxAF <- getMaxAF(hypothesis) 

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

    # If would return negative likelihood skip
    if(hypothesis$doDropin & checkDropin(x, maxAF))
	{
	if(logObjective) result = log10(0) else result = 0
	if(verbose) print(result)
	return(-result)
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
    -result
  }
  
  # Number of iterations
  if(is.null(iterMax)) iterMax = 1000
  
  lower = lower.bounds(args, zero, logDegradation)
  upper = upper.bounds(args, zero)
  lower = lower[names(template)] 
  upper = upper[names(template)] 


  list(#par     = unlist(template), 
       fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       #control = list(fnscale=-1, factr=1e7, maxit=500), 
       control = list(strategy=3,NP=4*length(unlist(upper)),itermax=iterMax) 
       #method  = "L-BFGS-B",
       #hessian = FALSE )
       )
}



getMaxAF <- function(hypothesis)
	{
	# Get the maximum allele frequency of any alleles in database 
	#
	# Parameters:
	# hypothesis: Hypothesis from which to get the maximum allele fraction
	maxAF <- function (db)
		{
		max(db[,1])
		}
	out <- max(sapply(X=hypothesis$alleleDb, FUN=maxAF))
	return(out)
	}


checkDropin <- function(params, maxAF)
	{
	# Check that dropin value will not create negative likelihoods
	#
	# Parameters:
	# 	params: parameters to be checked
	# 	maxAF: maximum allele fraction
	dropin = params$dropin
	dropout = params$dropout
	check = maxAF * (dropin*(1-dropout))
	if(any(check>1)) return(TRUE) else return(FALSE)
	}


DEoptimLoop <- function(PARAMS, tolerance=1e-6)
	{
	# Optimises over parameters, while checking for convergence every 50 iterations
	#
	# Parameters:
	# PARAMS: parameters for DEoptim generated by optimization.params
	# Set number of iterations to 50
	
	PARAMS$control$itermax = 50
	# Bogus result to check against at first
	oldresult = 999
	bestmemitOut <- NULL
	bestvalitOut <- NULL
	iterOut <- NULL
	nfevalOut <- NULL
	# Set check flag to FALSE
	flag = FALSE
	while(!flag)
		{
		# Run DEoptim for 100 iterations
		results <- do.call(DEoptim,PARAMS)
		# Get results for output
		bestmemitOut <- rbind(bestmemitOut, results$member$bestmemit)
		bestvalOut <- rbind(bestvalitOut, results$member$bestvalit)
		iterOut <- c(iterOut, results$optim$iter)
		nfevalOut <- c(results$optim$nfeval) 
		# Check if result is the same as the pervious one (50 iterations before)
		condition = abs((results$optim$bestval-oldresult)/oldresult)<tolerance
		if(condition) flag=TRUE
		# If not set the initial population to the last population from the last result
		PARAMS$control$initialpop = results$member$pop
		# Set the result to be checked against as the current result
		oldresult = results$optim$bestval
		}
	# Update some results to include all steps of loop
	results$member$bestmemit = bestmemitOut
	results$member$bestvalit = bestvalitOut
	results$optim$iter = sum(iterOut)
	results$optim$nfeval = sum(nfevalOut)
	return(results)
	}

