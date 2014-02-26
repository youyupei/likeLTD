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
                                withPenalties=TRUE, objective=NULL, iterMax=NULL, likeMatrix=FALSE,...) {
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
  if(is.null(objective)) objective = create.likelihood.vectors(hypothesis, likeMatrix=likeMatrix,...)

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
      x = relistArguments(x, hypothesis, fixed=fixed, arguments=arguments, logDegradation=logDegradation)
    # Otherwise, checks some options.
    else { 
      # Make sure it contains fixed terms.
      if(length(fixedstuff)) x = append(x, fixedstuff)
      # Converts degradation from log10 format.
      if(logDegradation && "degradation" %in% names(x))
        x$degradation = 10^x$degradation
    }

    # If would return negative likelihood skip
    if(hypothesis$doDropin & checkDropin(x, maxAF, hypothesis$nUnknowns+nrow(hypothesis$dropoutProfs)))
	{
	if(logObjective) result = log10(0) else result = 0
	if(verbose) print(result)
	return(-result)
	}

    # Compute objective function.
    result <- do.call(objective, x)
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
    # return result
    -result
  }
  
  # Number of iterations
  if(is.null(iterMax)) iterMax = 1000
  
  lower = lower.bounds(args, zero, logDegradation)
  upper = upper.bounds(args, zero)
  lower = lower[names(template)] 
  upper = upper[names(template)] 

 # population size for optimisation
  searchPop = 4*length(unlist(upper))
  # increases for relatedness
  searchPop = round(searchPop * relFactor(hypothesis$relatedness))


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


# factor by which to increase the population size for optimisation when taking into account relatedness
# optimisation seems more difficult when related, so need a more thorough search
relFactor = function(relatedness,base=1,fac1=9,fac2=4,fac3=2)
	{
	base + fac1*relatedness[1] + fac2*relatedness[2] + fac3*prod(relatedness)
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


checkDropin <- function(params, maxAF, nContrib)
	{
	# Check that dropin value will not create negative likelihoods
	#
	# Parameters:
	# 	params: parameters to be checked
	# 	maxAF: maximum allele fraction
	#	nContrib: number of unknowns plus number of contributors subject to dropout
	dropin = params$dropin
	dropout = params$dropout
	# dropin rate is just dropin probability if no reference individual
	if(nContrib==0)
		{
		check = maxAF * dropin
		} else {
		check = maxAF * (dropin*(1-dropout))
		}
	if(any(check>1)) return(TRUE) else return(FALSE)
	}

DEoptimLoop <- function(PARAMS, tolerance=1e-6)
	{
	# Optimises over parameters, while checking for convergence every 50 iterations
	#
	# Parameters:
	# 	PARAMS: parameters for DEoptim generated by optimization.params
	# 	Set number of iterations to 50
	
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



get.likely.genotypes = function(hypothesis,params,results,joint=FALSE,prob=0.1)
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
	locusCentricHyp = transform.to.locus.centric(hypothesis)
	# function to return genotype combinations for each locus
	genCombs = function(locusHyp,alleleDb)
		{
		genotypeIndex = likelihood.constructs.per.locus(locusHyp)$genotypes
		genotypes = matrix(rownames(alleleDb)[genotypeIndex],ncol=ncol(genotypeIndex))
		return(genotypes)
		}
	genotypes = mapply(FUN=genCombs,locusHyp=locusCentricHyp,alleleDb=hypothesis$alleleDb)
	# make a new output function to output genotype likelihoods
	newParams = optimisation.params(hypothesis,likeMatrix=TRUE)
	# get genotype likelihoods with parameters returned by previous optimisation (results)
	likes = newParams$fn(results$optim$bestmem)$objectives
	# convert to probabilities
	likes = lapply(likes,FUN=function(x) x/sum(x))
	# if we only want the joint distributions
	if(joint==TRUE) 
		{
		# joint genotypes
		outJoint = mapply(FUN=subGens,genotypes,likes,rotate=TRUE,SIMPLIFY=FALSE)
		# top genotype combination
		topGenotypes = mapply(FUN=function(a,b) a[,which.max(b)],genotypes,likes,SIMPLIFY=TRUE)
		return(list(genotypes=genotypes,probabilities=outJoint,topGenotypes=topGenotypes))
		}
	# function to get marginal probabilities and subset to those with prob>prob
	marginalProbs = function(genotypes,probabilities,indiv=1,prob=0.1,top=FALSE)
		{
		marginals = marginal(genotypes,probabilities,indiv)
		if(top==TRUE)	
			{
			topMarginals = marginals$genotypes[which.max(marginals$probabilities),,drop=FALSE]
			return(topMarginals)
			}
		subMarginals = subGens(marginals$genotypes,marginals$probabilities,prob)
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
	topGenotypes = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,genotypes=genotypes,probabilities=likes,indiv=x,prob=prob,top=TRUE,SIMPLIFY=TRUE),simplify=FALSE)
	topGenotypes = topGenotypes[index]	
	#return(list(genotypes=genotypes,probabilities=likes))
	return(list(marginals=out,topGenotypes=topGenotypes))
	}


# function to get marginal probabilities for a single contributor
# for internal use within getLikes()
marginal = function(genotypes,probabilities,indiv=1)
	{
	# genotypes for just 1 contributor
	allGen1 = genotypes[(((2*indiv)-1):(2*indiv)),]
	# remove duplicate genotypes
	outGens = unique(t(allGen1))
	# sum probabilities for each identical single genotype
	outProbs = apply(outGens,MARGIN=1,FUN=function(x) sum(probabilities[getMatching(singleGens=allGen1,matchGen=x)]))
	return(list(genotypes=outGens,probabilities=outProbs))
	}

# function to return an index of which genotypes in group match a given genotype
# for internal use in marginal()
getMatching = function(singleGens,matchGen)
	{
	index1 = colSums(apply(singleGens,MARGIN=2,FUN=function(x) x%in%matchGen))
	index2 = colSums(apply(singleGens,MARGIN=2,FUN=function(x) matchGen%in%x))
	outIndex = which((index1+index2)==4)
	return(outIndex)
	}

# function to subset genotypes to those with probability greater than prob
# for internal use within getLikes()
subGens = function(genotypes,probabilities,rotate=FALSE,prob=0.1)
	{
	if(rotate==TRUE) genotypes = t(genotypes)
	genotypes = genotypes[rev(order(probabilities)),]
	probabilities = sort(probabilities,decreasing=TRUE)
	index = which(probabilities>prob)
	genotypes = genotypes[index,]
	probabilities = probabilities[index]
	return(list(genotypes=genotypes,probabilities=probabilities))
	}


