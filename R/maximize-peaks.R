upper.bounds.peaks = function(arguments, nloc, zero=1e-6, logDegradation=FALSE) { 
  # Upper bounds of optimisation function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimisation function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  degradation = if(logDegradation) { 0-zero } else { 1-zero }
  degradation = rep(degradation, length(arguments$degradation))
  rcont       = rep(100, length(arguments$rcont))
  sigma        = 10000
  dropin      = NULL
  stutter     = rep(.08-zero,nloc)
  if(!is.null(arguments[["dropin"]])) dropin = 10 - zero

  list(degradation     = degradation,
       rcont           = rcont,
       sigma           = sigma,
       stutter         = stutter,
       dropin          = dropin)[names(arguments)]
}


lower.bounds.peaks = function(arguments, nloc, zero=1e-6, logDegradation=FALSE) { 
  # Lower bounds of optimisation function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimisation function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  #   logDegradation: Wether degradation parameters are entered as exponents of
  #                   10.
  degradation = if(logDegradation) { -20 } else { 0 }
  degradation = rep(degradation, length(arguments$degradation))
  rcont       = rep(zero, length(arguments$rcont))
  sigma       = 1
  stutter     = rep(.08-zero,nloc)
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = zero

  list(degradation     = degradation,
       rcont           = rcont, 
       sigma           = sigma,
       stutter         = stutter,
       dropin          = dropin)[names(arguments)]
}


optimisation.params.peaks <- function(hypothesis, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=1e-6, throwError=FALSE,
                                withPenalties=TRUE, objective=NULL, iterMax=75, likeMatrix=FALSE,...) {
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

  hypothesis = likeLTD:::add.args.to.hypothesis(hypothesis, ...)
  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  # If the objective function has not been handed to optimizatio.params,
  # make the objective function
  if(is.null(objective)) objective = create.likelihood.vectors.peaks(hypothesis, likeMatrix=likeMatrix,...)


  # Get maximum allele fraction
  maxAF <- likeLTD:::getMaxAF(hypothesis) 


  args = arguments
  if(is.null(args)) args = initial.arguments.peaks(hypothesis)
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
      x = relistArguments.peaks(x, hypothesis, fixed=fixed, arguments=arguments, logDegradation=logDegradation)
    # Otherwise, checks some options.
    else { 
      # Make sure it contains fixed terms.
      if(length(fixedstuff)) x = append(x, fixedstuff)
      # Converts degradation from log10 format.
      if(logDegradation && "degradation" %in% names(x))
        x$degradation = 10^x$degradation
    }

    # If would return negative likelihood skip
    if(hypothesis$doDropin)
        {
        if(likeLTD:::checkDropin(x, maxAF, hypothesis$nUnknowns+nrow(hypothesis$knownProfs)))
	        {
	        if(logObjective) result = log10(0) else result = 0
	        if(verbose) print(result)
	        return(-result)
	        }
	    }

    # Compute objective function.
    result <- do.call(objective, x)

    tmpResult = result

   #if(any(is.infinite(result$objectives))) TOTALRES <<- result

if(likeMatrix==TRUE) return(result)
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
    if(is.na(result))
        {
        OVERALLRES <<- tmpResult
        }
    # if result is infinite make sure it returns -Inf
    if(is.infinite(result)|is.na(result)) result = -Inf
    # return result
    -result
  }


  
  lower = lower.bounds.peaks(args, ncol(hypothesis$queriedProfile), zero, logDegradation)
  upper = upper.bounds.peaks(args, ncol(hypothesis$queriedProfile), zero, logDegradation)
  lower = lower[names(template)] 
  upper = upper[names(template)] 



 # population size for optimisation
  searchPop = 4*length(unlist(upper))


  list(#par     = unlist(template), 
       fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       #control = list(fnscale=-1, factr=1e7, maxit=500), 
       control = list(strategy=3,NP=searchPop,itermax=iterMax) 
       #method  = "L-BFGS-B",
       #hessian = FALSE )
       )
}


# Best(?) guess for initial arguments. 
initial.arguments.peaks <- function(hypothesis, ...) {
  # Best(?) guess for initial arguments. 
  #
  # Parameters: 
  #    hypothesis: Hypothesis for which to guess nuisance paramters. 

  hypothesis = likeLTD:::add.args.to.hypothesis(hypothesis, ...)
  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  # -1 because relative to first.
  nrcont          = max(nrow(hypothesis$knownProfs)
                        + hypothesis$nUnknowns - 1, 0)
  degradation     = rep( 3e-3, 
                         nrow(hypothesis$knownProfs) + hypothesis$nUnknowns )
  rcont           = runif(nrcont, min=0.5, max=1.5)
  dropin          = NULL
  sigma           = 1/4
  stutter         = rep(0.08,times=ncol(hypothesis$queriedProfile))
  if(hypothesis$doDropin) dropin = 1e-2


  list(degradation     = degradation,
       rcont           = rcont,
       sigma           = sigma,
       stutter         = stutter,
       dropin          = dropin)
}


relistArguments.peaks <- function( parameters, hypothesis, fixed=NULL,
                             logDegradation=TRUE, arguments=NULL ) {
  # Remakes arguments from a flat vector into a list.
  if(is.null(arguments)) arguments = initial.arguments.peaks(hypothesis) 
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


get.likely.genotypes.peaks = function(hypothesis,params,results,joint=FALSE,prob=ifelse(joint==FALSE,0.1,0.05))
	{
	# Function to return likely genotypes for each individual
	# Finds the marginal probabilities of genotypes, and then
	# returns the most likely genotypes, with their probabilities
	#
	# Parameters:
	#	hypothesis: hypothesis from defence.hypothesis(args)
	#	params: parameters object returned from optimisation.params(hypothesis)
	#	results: results object returned from DEoptimLoop(params)
	#	prob: genotypes with a probability greater than this will be returned

	# transform hypothesis to locus centric
	locusCentricHyp = transform.to.locus.centric.peaks(hypothesis)
	# function to return genotype combinations for each locus
	genCombs = function(locusHyp,alleleDb)
		{
		genotypeIndex = likelihood.constructs.per.locus.peaks(locusHyp)$genotypes
		genotypes = matrix(rownames(alleleDb)[genotypeIndex],ncol=ncol(genotypeIndex))
		return(genotypes)
		}
	genotypes = mapply(FUN=genCombs,locusHyp=locusCentricHyp,alleleDb=hypothesis$alleleDb)
	# make a new output function to output genotype likelihoods
	newParams = optimisation.params.peaks(hypothesis,likeMatrix=TRUE)
	# get genotype likelihoods with parameters returned by previous optimisation (results)
	likes = newParams$fn(results$optim$bestmem)$objectives
	# convert to probabilities
	likes = lapply(likes,FUN=function(x) x/sum(x))
	# if we only want the joint distributions
	if(joint==TRUE) 
		{
		# joint genotypes
		outJoint = mapply(FUN=likeLTD:::subGens,genotypes,likes,rotate=TRUE,prob=prob,SIMPLIFY=FALSE)
		# top genotype combination
		topGenotypes = mapply(FUN=function(a,b) a[,which.max(b)],genotypes,likes,SIMPLIFY=TRUE)
		topGenotypes = t(topGenotypes)
		# get top probability
		topProbability = prod(sapply(likes,FUN=max))
		return(list(joint=outJoint,topGenotypes=list(genotype=topGenotypes,probability=topProbability)))
		}
	# function to get marginal probabilities and subset to those with prob>prob
	marginalProbs = function(genotypes,probabilities,indiv=1,prob=0.1,top=FALSE)
		{
		marginals = likeLTD:::marginal(genotypes,probabilities,indiv)
		if(top==TRUE)	
			{
			topMarginals = marginals$genotypes[which.max(marginals$probabilities),,drop=FALSE]
			topProbability = max(marginals$probabilities)		
			return(list(genotype=topMarginals,probability=topProbability))
			}
		subMarginals = likeLTD:::subGens(marginals$genotypes,marginals$probabilities,prob)
		return(subMarginals)
		}
	# number of contributors
	ncont = nrow(genotypes[[1]])/2
	# get marginal and subset at every locus for every individual
	out = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,genotypes=genotypes,probabilities=likes,indiv=x,prob=prob,SIMPLIFY=FALSE),simplify=FALSE)
	# order by dropout rate
	rcont = vector(length=ncont)
	rcont[hypothesis$refIndiv] = 1
	rcont[-hypothesis$refIndiv] = results$optim$bestmem[grep("rcont",names(results$optim$bestmem))]
	index = (1:ncont)[rev(order(rcont))]
	out = out[index]
	# get top genotypes for marginals
	topGenotypes = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,genotypes=genotypes,probabilities=likes,indiv=x,prob=prob,top=TRUE,SIMPLIFY=FALSE),simplify=FALSE)
	# get top probabilities for marginals	
	topProbabilities = sapply(1:ncont,FUN=function(y) prod(sapply(topGenotypes[[y]],FUN=function(x) x$probability)),simplify=FALSE)
	topGenotypes = sapply(1:ncont,FUN=function(y) sapply(topGenotypes[[y]],FUN=function(x) x$genotype),simplify=FALSE)
	topGenotypes = topGenotypes[index]	
	topGenotypes = lapply(topGenotypes,FUN=t)
	topProbabilities = topProbabilities[index]
	#return(list(genotypes=genotypes,probabilities=likes))
	return(list(out,topGenotypes=list(genotypes=topGenotypes,probabilities=topProbabilities)))
	}

