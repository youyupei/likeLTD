upper.bounds.peaks = function(arguments, nloc, zero=1e-6, logDegradation=FALSE) { 
  # Upper bounds of optimisation function.
  # 
  # Parameters:
  #   arguments: Arguments passed to the optimisation function. Used as a
  #              template.
  #   zero: Some bounds should be given as >, rather than >=. This arguments is
  #         an epsilon to simulate the first case.
  degradation = if(logDegradation) { -1 } else { 1-zero }
  degradation = rep(degradation, length(arguments$degradation))
  DNAcont       = rep(5000, length(arguments$DNAcont))
  scale        = 1000
  gradientS = 0.01
  gradientAdjust     = rep(5,nloc)
  #interceptAdjust     = rep(5,nloc)
#locusAdjust     = rep(5,nloc)
  #interceptS = 1e-2
  meanD = NULL
  if(!is.null(arguments[["meanD"]])) meanD = 0.1
  meanO = NULL
  if(!is.null(arguments[["meanO"]])) meanO = 0.1
  repAdjust   = rep(5,length(arguments$repAdjust))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = 100

  list(degradation     = degradation,
       DNAcont           = DNAcont,
       scale           = scale,
       gradientS = gradientS,
       gradientAdjust         = gradientAdjust,
#locusAdjust     = locusAdjust,
       #interceptAdjust         = interceptAdjust,
       #interceptS   = interceptS,
       repAdjust       = repAdjust,
       meanD = meanD,
       meanO = meanO,
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
  DNAcont       = rep(zero, length(arguments$DNAcont))
  scale       = 0+zero
  gradientS = 0+zero
  gradientAdjust     = rep(0.2,nloc)
  #interceptAdjust     = rep(0.2,nloc)
#locusAdjust     = rep(0.2,nloc)
  #interceptS  = 1e-15
  meanD = NULL
  if(!is.null(arguments[["meanD"]])) meanD = 0+zero
  meanO = NULL
  if(!is.null(arguments[["meanO"]])) meanO = 0+zero
  repAdjust   = rep(0.2,length(arguments$repAdjust))
  dropin      = NULL
  if(!is.null(arguments[["dropin"]])) dropin = 5

  list(degradation     = degradation,
       DNAcont           = DNAcont, 
       scale           = scale,
       gradientS = gradientS,
       gradientAdjust         = gradientAdjust,
      # interceptAdjust         = interceptAdjust,
#locusAdjust         = locusAdjust,
      # interceptS = interceptS,
       repAdjust       = repAdjust,
       meanD = meanD,
       meanO = meanO,
       dropin          = dropin)[names(arguments)]
}


optimisation.params.peaks <- function(hypothesis, verbose=TRUE, fixed=NULL,
                                logObjective=TRUE, logDegradation=TRUE,
                                arguments=NULL, zero=1e-6, throwError=FALSE,
                                withPenalties=TRUE, doLinkage=TRUE, objective=NULL, iterMax=75, 
				likeMatrix=FALSE,diagnose=FALSE,...) {
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
  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  # If the objective function has not been handed to optimizatio.params,
  # make the objective function
  if(is.null(objective)) objective = create.likelihood.vectors.peaks(hypothesis, likeMatrix=likeMatrix, diagnose=diagnose,...)
  # Get maximum allele fraction
  maxAF <- getMaxAF(hypothesis) 
  args = arguments
  if(is.null(args)) args = initial.arguments.peaks(hypothesis)
  if(logDegradation && "degradation" %in% names(args)) 
    args$degradation = log10(args$degradation)
  # Make sure we don't include empty stuff (like rcont sometimes)
  template = Filter(function(n) length(n) > 0, args)
  # deal with fixed variables
  if(!is.null(fixed)) {
    fixedstuff = args[fixed]
    template = args[-which(names(args) %in% fixed)]
  } else  fixedstuff = NULL
    # Linkage adjustment if brothers
    linkBool = doLinkage&identical(hypothesis$relatedness,c(0.5,0.25))&hypothesis$hypothesis=="prosecution"
  if(linkBool)
    {
    linkFactor = linkage(hypothesis)
    if(logObjective) linkFactor = log10(linkFactor)
    }
    # result function
    result.function <- function(x) 
        {
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
        if(checkDropin(x, maxAF, hypothesis$nUnknowns+nrow(hypothesis$knownProfs)))
	        {
	        if(logObjective) result = log10(0) else result = 0
	        if(verbose) print(result)
	        return(-result)
	        }
	    }
	# if stutter is >100% or <0% return negative likelihood
    stuttermodel = function(stutterMean,stutterGradient,i)
    	{
    	stutterMean+(stutterGradient*i)
    	}
    toAdd = 0
    if(hypothesis$doDoubleStutter) toAdd = toAdd + x$meanD
    if(hypothesis$doOverStutter) toAdd = toAdd + x$meanO
#    condition1 = mapply(x$gradientAdjust*x$gradientS,0,hypothesis$alleleDb,
#                FUN=function(gradA,stuttA,db) 
#                    any(stuttermodel(stuttA,gradA,db[,3])+toAdd>1)|any(stuttermodel(stuttA,gradA,db[,3])+toAdd<0))
#    condition2 = mapply(x$gradientAdjust*x$gradientS,0,hypothesis$alleleDb,
#                FUN=function(gradA,stuttA,db) 
#                    any(stuttermodel(stuttA,gradA,db[,3])>1)|any(stuttermodel(stuttA,gradA,db[,3])<0))
condition1 = mapply(x$gradientAdjust*x$gradientS,hypothesis$alleleDb,
                FUN=function(gradA,db) 
                    any(stuttermodel(0,gradA,db[,3])+toAdd>1)|any(stuttermodel(0,gradA,db[,3])+toAdd<0))
    condition2 = mapply(x$gradientAdjust*x$gradientS,hypothesis$alleleDb,
                FUN=function(gradA,db) 
                    any(stuttermodel(0,gradA,db[,3])>1)|any(stuttermodel(0,gradA,db[,3])<0))
    if(any(condition1)|any(condition2))
	    {
	    if(logObjective) result = log10(0) else result = 0
	    if(verbose) print(result)
	    return(-result)
	    }
    # Compute objective function.
    result <- do.call(objective, x)

	if(likeMatrix==TRUE|diagnose==TRUE) return(result)

    # Compute as log if requested, otherwise as product.
	if(withPenalties) 
		{
		if(logObjective)
			{
			result <- sum(log10(result$objectives) + log10(result$penalties))
			} else {
			result <- prod(result$objectives) * prod(result$penalties)
			}
		} else {
		if(logObjective) {
			result <- sum(log10(result$objectives))
			} else {
			result <- prod(result$objectives)
			}
		} 
		# Print out if requested.
		if(verbose) 
			{
			# print(unlist(append(x, list(result=result))))
			print(result)
			}
		# If result is infinite, do throw
		if(throwError == TRUE && is.infinite(result)) 
			{
			cat("Objective function is over/underflow: ", result, "\n")
			print(x)
			stop("Objective function over/underflow")
			}
		# if result is infinite make sure it returns -Inf
		if(is.infinite(result)|is.na(result)) result = -Inf

    		# Linkage adjustment
		if(linkBool)
		    {
		    if(logObjective==TRUE) 
		        {
		        result = result+linkFactor
		        } else {
		        result = result*linkFactor
		        }
		    }
		# return result
		-result
	}
  lower = lower.bounds.peaks(args, ncol(hypothesis$queriedProfile), zero, logDegradation)
  upper = upper.bounds.peaks(args, ncol(hypothesis$queriedProfile), zero, logDegradation)
  lower = lower[names(template)] 
  upper = upper[names(template)] 
 # population size for optimisation
  searchPop = 4*length(unlist(upper))
  list(fn      = result.function, 
       lower   = unlist(lower), 
       upper   = unlist(upper),
       control = list(strategy=3,NP=searchPop,itermax=iterMax) 
       )
}


# Best(?) guess for initial arguments. 
initial.arguments.peaks <- function(hypothesis, ...) {
  # Best(?) guess for initial arguments. 
  #
  # Parameters: 
  #    hypothesis: Hypothesis for which to guess nuisance paramters. 

  hypothesis = add.args.to.hypothesis(hypothesis, ...)
  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  # -1 because relative to first.
  nDNAcont          = max(nrow(hypothesis$knownProfs)
                        + hypothesis$nUnknowns, 0)
  degradation     = rep( 3e-3, 
                         nrow(hypothesis$knownProfs) + hypothesis$nUnknowns )
  DNAcont           = runif(nDNAcont, min=0.5, max=1.5)
  dropin          = NULL
  meanD    = NULL
  meanO    = NULL
  scale           = 1/4
  gradientS = 0.08
  gradientAdjust   = rep(1,times=ncol(hypothesis$queriedProfile))
 # interceptAdjust   = rep(1,times=ncol(hypothesis$queriedProfile))
#locusAdjust   = rep(1,times=ncol(hypothesis$queriedProfile))
 # interceptS = 0
  repAdjust       = rep(1,times=max(length(hypothesis$peaksProfile)-1,0))
  if(hypothesis$doDropin) dropin = 20
  if(hypothesis$doDoubleStutter) meanD = 0.02
  if(hypothesis$doOverStutter) meanO = 0.02


  list(degradation     = degradation,
       DNAcont           = DNAcont,
       scale           = scale,
       gradientS = gradientS,
       gradientAdjust         = gradientAdjust,
   #    interceptAdjust         = interceptAdjust,
#       locusAdjust         = locusAdjust,
   #    interceptS = interceptS,
       repAdjust       = repAdjust,
       meanD = meanD,
       meanO = meanO,
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


get.likely.genotypes.peaks = function(hypothesis,params,results,posterior=FALSE,joint=FALSE,prob=ifelse(joint==FALSE,0.1,0.05))
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
	locusCentricHyp = likeLTD:::transform.to.locus.centric.peaks(hypothesis)
	# function to return genotype combinations for each locus
	genCombs = function(locusHyp,alleleDb)
		{
		genotypes = likeLTD:::likelihood.constructs.per.locus.peaks(locusHyp)$genotypes
		return(genotypes)
		}
	genotypes = mapply(FUN=genCombs,locusHyp=locusCentricHyp,alleleDb=hypothesis$alleleDb,SIMPLIFY=FALSE)
	# make a new output function to output genotype likelihoods
	newParams = optimisation.params.peaks(hypothesis,likeMatrix=TRUE)
	# get genotype likelihoods with parameters returned by previous optimisation (results)
	sepLikes = newParams$fn(results$optim$bestmem)$objectives
	combinedLikes = sapply(1:ncol(sepLikes),FUN=function(x) sepLikes[1,x][[1]]*sepLikes[2,x][[1]],simplify=FALSE)
	names(combinedLikes) = colnames(sepLikes)
	# convert to probabilities
	likes = lapply(combinedLikes,FUN=function(x) x/sum(x))
	# if want posterior, return
	if(posterior==TRUE) return(list(genotypes=genotypes,probabilities=likes))
	# if we only want the joint distributions
	if(joint==TRUE) 
		{
		# joint genotypes
		outJoint = mapply(FUN=subGens,genotypes,likes,rotate=TRUE,prob=prob,SIMPLIFY=FALSE)
		# top genotype combination
		topGenotypes = mapply(FUN=function(a,b) a[,which.max(b)],genotypes,likes,SIMPLIFY=TRUE)
		topGenotypes = t(topGenotypes)
		# get top probability
		topProbability = prod(sapply(likes,FUN=max))
		return(list(joint=outJoint,topGenotypes=list(genotype=topGenotypes,probability=topProbability)))
		}
	# function to get marginal probabilities and subset to those with prob>prob
	marginalProbs = function(gens,probs,indiv=1,prob=0.1,top=FALSE)
		{
		marginals = marginal(gens,probs,indiv)
		if(top==TRUE)	
			{
			topMarginals = marginals$genotypes[which.max(marginals$probabilities),,drop=FALSE]
			topProbability = max(marginals$probabilities)		
			return(list(genotype=topMarginals,probability=topProbability))
			}
		subMarginals = subGens(marginals$genotypes,marginals$probabilities,prob)
		return(subMarginals)
		}
	# number of contributors
	ncont = nrow(genotypes[[1]])/2
	# get marginal and subset at every locus for every individual
	out = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,gens=genotypes,probs=likes,indiv=x,prob=prob,SIMPLIFY=FALSE),simplify=FALSE)
	# get top genotypes for marginals
	topGenotypes = sapply(1:ncont,FUN=function(x) mapply(FUN=marginalProbs,gens=genotypes,probs=likes,indiv=x,prob=prob,top=TRUE,SIMPLIFY=FALSE),simplify=FALSE)
	# get top probabilities for marginals	
	topProbabilities = sapply(1:ncont,FUN=function(y) prod(sapply(topGenotypes[[y]],FUN=function(x) x$probability)),simplify=FALSE)
	topGenotypes = sapply(1:ncont,FUN=function(y) sapply(topGenotypes[[y]],FUN=function(x) x$genotype),simplify=FALSE)
	topGenotypes = lapply(topGenotypes,FUN=t)
	#return(list(genotypes=genotypes,probabilities=likes))
	return(list(out,topGenotypes=list(genotypes=topGenotypes,probabilities=topProbabilities)))
	}


plot.peaks.results = function(hyp,res,replicate=1,fileName="peakHeights.pdf")
	{
	# mean & sd from results
	diagParams = optimisation.params.peaks(hyp,diagnose=TRUE)
	info = diagParams$fn(res$optim$bestmem)
	# genotype likelihoods from results
	likesParams = optimisation.params.peaks(hyp,likeMatrix=TRUE)
	likes = likesParams$fn(res$optim$bestmem)
	maxIndex = sapply(likes$objectives, FUN=which.max)

	# plot for each locus
	pdf(fileName)
	par(mfrow=rep(ceiling(sqrt(length(hyp$alleleDb))),times=2),mar=c(3,2,2,0))
	for(i in 1:length(hyp$alleleDb))
		{
		# results from top genotype at this locus
		if(is.matrix(info[[i]])) 
			{
			topGenoInfo = info[[i]][maxIndex[i],replicate][[1]]
			} else {
			topGenoInfo = info[[i]]
			}
		if(length(topGenoInfo)==1) topGenoInfo = topGenoInfo[[1]]
		
		heights = topGenoInfo$height[order(as.numeric(names(topGenoInfo$height)))]
		scale = res$optim$bestmem[grep("scale",names(res$optim$bestmem))]
		shapes = topGenoInfo$mu/scale
		shapes = shapes[order(as.numeric(names(shapes)))]
		#  confidence intervals
		CIs = sapply(shapes,FUN=function(x) qgamma(p=c(0.025,0.25,0.5,0.75,0.975),shape=x,scale=scale))
		YLIM = c(0,max(c(CIs,heights)))
		# plot
		boxplot(CIs,ylim=YLIM,range=0,main=names(hyp$alleleDb)[i])
		boxplot(t(heights),ylim=YLIM,border="red",add=TRUE)
		}
	dev.off()

	}


# check for convergence
checkConverged = function(new,old,tol)
    {
    return(abs((new-old)/old)<tol)
    }
multiConverged = function(L,globalVal,tolerance)
    {
    any(sapply(L[(length(L)-5):(length(L)-1)],FUN=function(x) !checkConverged(L[length(L)],x,tolerance)))|!checkConverged(L[length(L)],globalVal,tolerance)
    }

evaluate.peaks <- function(P.pars, D.pars, tolerance=1e-6, n.steps=NULL, scaleLimit=1, interim=FALSE, CR.start=0.1, CR.end=0.7, seed.input=NULL){
	# P.pars D.pars: parameter object created by optimisation.params()
	# the smallest convergence threshold (ie for the last step)
	# number of convergence thresholds
	
	# for each step, run a DEoptimLoop both for P and D, until each converges at that steps accuracy
	if(is.null(seed.input)) 
	    {
	    seed.used =  as.numeric(Sys.time())
	    } else {
        seed.used = seed.input
	    }
    set.seed(seed.used)

	# combine the outputs outside the loop
	P.bestmemit <- D.bestmemit <- NULL
	P.bestvalit <- D.bestvalit <- NULL
	P.iter <- D.iter <- NULL
	P.nfeval <- D.nfeval <- NULL

	# run first step
	n = 1

	# change DEoptim parameters
	P.pars$control$CR <- CR.start
	D.pars$control$CR <- CR.start
		
	# run DEoptimLoop until convergence at the required step
	D.step <- DEoptimLoop(D.pars,10)
	P.step <- DEoptimLoop(P.pars,10)

	# set global results
	GlobalDval = D.step$optim$bestval
	GlobalDmem = D.step$optim$bestmem
	GlobalPval = P.step$optim$bestval
	GlobalPmem = P.step$optim$bestmem

	# put the Def outputs into the combined output	
	D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
	D.bestvalit <- rbind(D.bestvalit, D.step$member$bestvalit)
	D.iter <- c(D.iter, D.step$optim$iter)
	D.nfeval <- c(D.nfeval, D.step$optim$nfeval)	

	# put the Pros outputs into the combined output	
	P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
	P.bestvalit <- rbind(P.bestvalit, P.step$member$bestvalit)
	P.iter <- c(P.iter, P.step$optim$iter)
	P.nfeval <- c(P.nfeval, P.step$optim$nfeval)	

	# recycle the current pop into the next loop
	D.pars$control$initialpop <- D.step$member$pop
	P.pars$control$initialpop <- P.step$member$pop

	# get standard mean standard deviation of initial optimisation phase
	#sdStep = mean(c(sd(P.step$member$bestvalit[1:75][!is.infinite(P.step$member$bestvalit[1:75])]),sd(D.step$member$bestvalit[1:75][!is.infinite(D.step$member$bestvalit[1:75])])))
	sdStep = max(c(sd(P.step$member$bestvalit[1:75][!is.infinite(P.step$member$bestvalit[1:75])]),sd(D.step$member$bestvalit[1:75][!is.infinite(D.step$member$bestvalit[1:75])])))
	# sometimes sd is very low (below 1 e.g. 3locus test)
	# if so set sd to >1 so log2(sd) is positive
	#if(sdStep<1) sdStep = 1.5
	# decide how many steps to run
	#if(is.null(n.steps)) n.steps = ceiling(log2(sdStep))*8+length(grep("cont",names(D.pars$upper)))
	if(is.null(n.steps)) n.steps = ceiling(sdStep+length(grep("cont",names(D.pars$upper))))

	# retain all the likelihood ratios
	Ld <- numeric(n.steps)
	Lp <- numeric(n.steps)

	Ld[n] = D.step$optim$bestval
	Lp[n] = P.step$optim$bestval

	# if more than one step, start loop
	if(n.steps>1)
		{
		# adjust tolerance gradually
		tol.steps <- geometric.series(10,tolerance,n.steps)
	
		# adjust DEoptim parameters gradually, so search space is confined more at the end
		CR.steps <- geometric.series(CR.start,CR.end,n.steps)

		for(n in 2:n.steps){
			# change DEoptim parameters
			D.pars$control$CR <- CR.steps[n]
			P.pars$control$CR <- CR.steps[n]
		
			# run DEoptimLoop until convergence at the required step
			D.step <- DEoptimLoop(D.pars,tol.steps[n])
            Ld[n] = D.step$optim$bestval
            P.step <- DEoptimLoop(P.pars,tol.steps[n])
            Lp[n] = P.step$optim$bestval

			# set global results
			if(D.step$optim$bestval<GlobalDval)
				{
				GlobalDval = D.step$optim$bestval
				GlobalDmem = D.step$optim$bestmem
				}
			if(P.step$optim$bestval<GlobalPval)
				{
				GlobalPval = P.step$optim$bestval
				GlobalPmem = P.step$optim$bestmem
				}

			# put the Def outputs into the combined output	
			D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
			D.bestvalit <- c(D.bestvalit, D.step$member$bestvalit)
			D.iter <- c(D.iter, D.step$optim$iter)
			D.nfeval <- c(D.nfeval, D.step$optim$nfeval)

			# put the Pros outputs into the combined output	
			P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
			P.bestvalit <- c(P.bestvalit, P.step$member$bestvalit)
			P.iter <- c(P.iter, P.step$optim$iter)
			P.nfeval <- c(P.nfeval, P.step$optim$nfeval)	

			# recycle the current pop into the next loop
			D.pars$control$initialpop <- D.step$member$pop
			P.pars$control$initialpop <- P.step$member$pop

			# generate outputs if interim = TRUE
			if(interim==TRUE) interim(P.step,D.step,n,n.steps)
		    }
		}
        # check for convergence   
        while(multiConverged(Ld,GlobalDval,tolerance)|multiConverged(Lp,GlobalPval,tolerance))
            {
	    n = n+1
            D.step <- DEoptimLoop(D.pars,tolerance)
            P.step <- DEoptimLoop(P.pars,tolerance)
            Ld = c(Ld,D.step$optim$bestval)
            Lp = c(Lp,P.step$optim$bestval)
			# set global results
			if(D.step$optim$bestval<GlobalDval)
				{
				GlobalDval = D.step$optim$bestval
				GlobalDmem = D.step$optim$bestmem
				}
			if(P.step$optim$bestval<GlobalPval)
				{
				GlobalPval = P.step$optim$bestval
				GlobalPmem = P.step$optim$bestmem
				}
			# put the Def outputs into the combined output	
			D.bestmemit <- rbind(D.bestmemit, D.step$member$bestmemit)
			D.bestvalit <- c(D.bestvalit, D.step$member$bestvalit)
			D.iter <- c(D.iter, D.step$optim$iter)
			D.nfeval <- c(D.nfeval, D.step$optim$nfeval)	
			# put the Pros outputs into the combined output	
			P.bestmemit <- rbind(P.bestmemit, P.step$member$bestmemit)
			P.bestvalit <- c(P.bestvalit, P.step$member$bestvalit)
			P.iter <- c(P.iter, P.step$optim$iter)
			P.nfeval <- c(P.nfeval, P.step$optim$nfeval)
			# recycle the current pop into the next loop
			D.pars$control$initialpop <- D.step$member$pop
			P.pars$control$initialpop <- P.step$member$pop
			# generate outputs if interim = TRUE
			if(interim==TRUE) interim(P.step,D.step,n,n.steps)
            }

    # get WoE
    WoE = Ld-Lp

	changeFlag=FALSE
	# if global result is better than result from last chunk
	if(P.step$optim$bestval>GlobalPval)
		{
		P.step$optim$bestval = GlobalPval
		P.step$optim$bestmem = GlobalPmem
		print("*** Final prosecution result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}
	if(D.step$optim$bestval>GlobalDval)
		{
		D.step$optim$bestval = GlobalDval
		D.step$optim$bestmem = GlobalDmem
		print("*** Final defence result was not the global optimum - consider re-running optimisation ***")
		changeFlag=TRUE
		}

	if(changeFlag) WoE <- c(WoE,D.step$optim$bestval - P.step$optim$bestval)


	# update final Pros results
	P.results <- P.step
	P.results$optim$bestval = GlobalPval
	P.results$optim$bestmem = GlobalPmem
	P.results$member$bestmemit <- P.bestmemit
	P.results$member$bestvalit <- P.bestvalit
	P.results$optim$iter <- P.iter
	P.results$optim$nfeval <- P.nfeval

	# update final Pros results
	D.results <- D.step
	D.results$optim$bestval = GlobalDval
	D.results$optim$bestmem = GlobalDmem
	D.results$member$bestmemit <- D.bestmemit
	D.results$member$bestvalit <- D.bestvalit
	D.results$optim$iter <- D.iter
	D.results$optim$nfeval <- D.nfeval

	print(Ld)
	print(Lp)

# return all results
return(list(Pros=P.results,Def=D.results, WoE=WoE, Lp=Lp, Ld=Ld, seed=seed))}




