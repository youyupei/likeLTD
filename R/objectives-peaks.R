


# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.vectors.peaks <- function(hypothesis, addAttr=FALSE, likeMatrix=FALSE, diagnose=FALSE, ...) {

  hypothesis = add.args.to.hypothesis(hypothesis, ...)

  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  locusCentric = transform.to.locus.centric.peaks(hypothesis)

  functions <- mapply(create.likelihood.per.locus.peaks, locusCentric,
                      MoreArgs=list(addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose))

create.likelihood.per.locus.peaks(locusCentric[[1]],addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose)


  #if(is.null(hypothesis$dropinPenalty)) hypothesis$dropinPenalty = 2 
  if(is.null(hypothesis$degradationPenalty)) hypothesis$degradationPenalty = 50 
  if(is.null(hypothesis$stutterPenalty)) hypothesis$stutterPenalty = 0.2


  likelihood.vectors <- function(degradation=NULL, DNAcont=NULL, scale=NULL, gradientS = NULL, gradientAdjust=NULL,stutterAdjust=NULL, dropin=NULL, meanS=NULL, 
                                meanD = NULL,meanO = NULL,
                                 repAdjust=NULL, detectionThresh=hypothesis$detectionThresh, degradationPenalty=hypothesis$degradationPenalty,stutterPenalty=hypothesis$stutterPenalty, ...) {
    # Call each and every function in the array.
    arguments = list(degradation=degradation, DNAcont=DNAcont,
                     scale = scale, repAdjust=repAdjust,
gradientS = gradientS,meanS=meanS, 
meanD = meanD,meanO=meanO,
detectionThresh = detectionThresh,
                     degradationPenalty=degradationPenalty, stutterPenalty=stutterPenalty, dropin=dropin)
    callme <- function(objective,grad,stut) {
    args = append(arguments, list(gradientAdjust=grad,stutterAdjust=stut))
      do.call(objective, args)
    }
    if(length(gradientAdjust) == 1) gradientAdjust = rep(gradientAdjust, length(functions))
    if(length(stutterAdjust) == 1) stutterAdjust = rep(stutterAdjust, length(functions))
#    if(setequal(names(stutter), colnames(hypothesis$cspProfile)))
#      stutterAdjust <- stutterAdjust[colnames(hypothesis$cspProfile)]
    objectives = mapply(callme, functions, gradientAdjust,stutterAdjust)
    arguments = append(arguments, list(...))
if(diagnose==TRUE) return(objectives)
    # calculate penalties
    pens <- do.call(penalties.peaks, append(arguments,list(gradientAdjust=gradientAdjust,stutterAdjust=stutterAdjust,nloc=ncol(hypothesis$queriedProfile))))
    list(objectives=objectives, penalties=pens)
  }
  if(addAttr) {
    attr(likelihood.vectors, "hypothesis") <- hypothesis
    attr(likelihood.vectors, "functions") <- functions
  }
  return(likelihood.vectors)
}


create.likelihood.per.locus.peaks <- function(hypothesis, addAttr=FALSE, likeMatrix = FALSE, diagnose=FALSE) {
  # Creates a likelihood function for a given hypothesis and locus
  #
  # A hypothesis is given by the number of unknown contributors, whether to model
  # dropin, so on and so forth.

  cons = likelihood.constructs.per.locus.peaks(hypothesis)
  doR = !is.null(hypothesis$doR) && hypothesis$doR == TRUE

  result.function <- function(scale,gradientS,gradientAdjust,stutterAdjust,
meanS,
meanD=NULL,meanO=NULL,
repAdjust=NULL,
                              degradation=NULL, DNAcont=NULL, 
			      detectionThresh = NULL, dropin=NULL, ...) {
    # Likelihood function for a given hypothesis and locus
    #
    # This function is specific to the hypothesis for which it was created.
    # It hides everything except the nuisance parameters over which to
    # optimize.
    #
    # Parameters:
    #   power: a scalar floating point value.
    #   dropout: the dropout rate for each replicate.
    #   degradation: relative degradation from each profiled individual in this
    #                hypothesis
    #   rcont: relative contribution from each profiled individual and unknown
    #          contributor in this hypothesis If it is one less than that, then
    #          an additional 1 is inserted at a position given by
    #          hypothesis$refIndiv. In general, this means either the queried
    #          individual (if subject to dropout and prosecution hypothesis) or
    #          the first individual subject to droput is the reference individual.
    #   ...: Any other parameter, e.g. for penalty functions. genotypeArray=cons$genotypes,
    #        These parameters are ignored here.
    # Returns: A scalar value giving the likelihood for this locus and
    #          hypothesis
    repAdjust = c(1,repAdjust)
    if(is.null(degradation)) degradation = c()
    if(length(DNAcont) != hypothesis$nUnknowns + ncol(cons$knownPresence))
      stop(sprintf("DNAcont should be %d long.",
                   hypothesis$nUnknowns + ncol(cons$knownPresence)))
    if(length(degradation) != hypothesis$nUnknowns +
                              ncol(cons$knownPresence))
      stop(sprintf("degradation should be %d long.",
                   hypothesis$nUnknowns + ncol(cons$knownPresence)))
    if(any(DNAcont < 0)) stop("found negative DNA contribution.")
    if(any(degradation < 0)) stop("found negative degradation parameter.")
    if(hypothesis$doDropin && is.null(dropin)) 
      stop("Model requires missing argument 'dropin'")
    else if(is.null(dropin)) dropin = 0
    if(hypothesis$doDropin && dropin < 0) 
      stop("Dropin rate should be between 0 and 1 (included).")

    if(diagnose==TRUE)
	{
	repRes <- peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
				scale=scale,gradientS = gradientS, gradientAdjust=gradientAdjust,stutterAdjust=stutterAdjust, meanS=meanS,
				meanD = meanD, meanO=meanO,
degradation=degradation, 
				repAdjust=repAdjust,detectionThresh=detectionThresh,doR=doR,diagnose=diagnose)
	return(repRes)
	}

    repRes <- matrix(peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
				scale=scale,gradientS = gradientS,gradientAdjust=gradientAdjust, stutterAdjust=stutterAdjust,meanS=meanS, 
				meanD = meanD,meanO=meanO,
degradation=degradation, 
				repAdjust=repAdjust,detectionThresh=detectionThresh,doR=doR),ncol=length(hypothesis$peaksProfile))

    # combine replicates
    if(ncol(repRes)>1)
	{
	#res = likeLTD:::prod.matrix.col(t(repRes))
	res = prod.matrix.row(repRes)
	} else {
	res = repRes	
	}

    # multiply by allele probs
    factorsRes = res*cons$factors




    # WHAT TO DO HERE FOR LOGS?
    # Figure out likelihood for good and return.
    if(likeMatrix==FALSE)
	{
    	return(sum(factorsRes))
	} else {
    	return(factorsRes)
	}
  }

  if(addAttr) {
    attr(result.function, "hypothesis") <- hypothesis
    attr(result.function, "constructs") <- cons
  }
  result.function
}


likelihood.constructs.per.locus.peaks = function(hypothesis) {
  # Creates the locus-specific data needed by each likehood function.
  #
  # Parameters:
  #   hypothesis: A hypothesis, for instance one returned by
  #               prosecution.hypothesis(...) or defence.hypothesis(...)
  alleles = rownames(hypothesis$alleleDb)
  if(is.null(alleles)) stop("Could not figure out alleles names.")
  alleles.vector = function(n) alleles %in% unlist(n)
  cspPresence     = apply(hypothesis$binaryProfile, 1, alleles.vector)
  knownPresence = apply(hypothesis$knownProfs, 1, alleles.vector)
  if(!is.matrix(cspPresence))
    cspPresence = matrix(ncol=0, nrow=length(alleles))
  if(!is.matrix(knownPresence))
    knownPresence = matrix(ncol=0, nrow=length(alleles))

  missingReps = apply(hypothesis$binaryProfile, 1, is.na)

#  genotypes <- compatible.genotypes.peaks(cspPresence, knownPresence, hypothesis$knownProfs,alleles,
 #                                   hypothesis$nUnknowns, hypothesis$doDropin,
  #                                  missingReps)

genotypes = explain.all.peaks(cspPresence,knownPresence,hypothesis$knownProfs,alleles,hypothesis$nUnknowns,hypothesis$peaksProfile,hypothesis$heightsProfile,hypothesis$doDoubleStutter,hypothesis$doOverStutter)

# get index of which alleles are from known contributors - do not want population allele probabilities for known contributors
if(nrow(hypothesis$knownProfs)>0) 
	{
	kIndex = (hypothesis$nUnknowns*2)+(1:(2*nrow(hypothesis$knownProfs)))
	} else {
	kIndex = c()
	}

# multiply by matching factor for known genotypes (1+fst, 1+2fst, ...)
if(length(kIndex>0)) {

  # Only take into account relatedness between Q and X under Hd 
  if(hypothesis$hypothesis=="defence")
	{
	# exclude known genotypes
  	factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             hypothesis$relatedness) 
	} else {
  	factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             c(0,0))
	}
 #factors = factors * prod(1+(kIndex*hypothesis$fst))
} else {
  # Only take into account relatedness between Q and X under Hd 
  if(hypothesis$hypothesis=="defence")
	{
	# exclude known genotypes
  	factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             hypothesis$relatedness) 
	} else {
  	factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             hypothesis$nUnknowns, hypothesis$doDropin,
                             hypothesis$queriedProfile,
                             c(0,0))
	}
}

genotypes = matrix(as.numeric(rownames(hypothesis$alleleDb))[genotypes],ncol=ncol(genotypes))

  list(cspPresence=cspPresence, knownPresence=knownPresence,
        missingReps=missingReps,
       genotypes=genotypes, factors=factors,
       freqMat=hypothesis$alleleDb[, 1])
}


# function to be called at each iteration of maximisation
peaks.probabilities = function(hypothesis,cons,
        DNAcont,scale,gradientS,gradientAdjust,
        stutterAdjust,meanS,meanD=NULL,meanO=NULL,
        degradation,repAdjust,detectionThresh,
        doR=FALSE,diagnose=FALSE)
    {
    if(hypothesis$doDropin==TRUE)
        {
        # no dropin currently
        stop("Cannot handle dropin")
        } else {
        # combine mean and adjustment
        locusGradient = gradientS*gradientAdjust
        locusStutter = meanS*stutterAdjust
	    if(doR==TRUE|diagnose==TRUE)
		    {
        	# probabilities for each replicate
        	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
        	        peak.heights.per.locus(
        	            genotypeArray=cons$genotypes,
						alleles=hypothesis$peaksProfile[[x]],
						heights=hypothesis$heightsProfile[[x]],
						DNAcont=DNAcont,
						gradientS = locusGradient,
						meanD=meanD,
						meanO=meanO,
						meanS=locusStutter,
						scale=scale,degradation=degradation,
						fragLengths=hypothesis$alleleDb[,2],
						LUSvals=hypothesis$alleleDb[,3],
						repAdjust=repAdjust[x],
						detectionThresh=detectionThresh,
						diagnose=diagnose))
		} else {
		if(!is.null(meanD)&!is.null(meanO))
			{
			# single, double and over stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
		    	        .Call(.cpp.probabilityPeaksSDO,
		    	            genotypeArray=cons$genotypes,
		                    alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                    heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                    DNAcont=rep(DNAcont,each=2), gradientS=locusGradient,
		                    meanD=meanD,meanO=meanO,
		                    meanS=locusStutter,
		                    scale=scale,degradation=rep(1+degradation,each=2),
		                    fragLengths=hypothesis$alleleDb[,2],
		                    fragNames=as.numeric(rownames(hypothesis$alleleDb)),
                            stutterIndex = hypothesis$alleleDb[,3],
		                    repAdjust=repAdjust[x],detectionThresh=detectionThresh))

			} else if(is.null(meanD)&is.null(meanO)) {
		   	# single stutter only
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
		    	        .Call(.cpp.probabilityPeaksS,
		    	            genotypeArray=cons$genotypes,
		                    alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                    heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                    DNAcont=rep(DNAcont,each=2),gradientS=locusGradient,
		                    meanS=locusStutter,
		                    scale=scale,degradation=rep(1+degradation,each=2),
		                    fragLengths=hypothesis$alleleDb[,2],
		                    fragNames=as.numeric(rownames(hypothesis$alleleDb)),
						    stutterIndex = hypothesis$alleleDb[,3],
		                    repAdjust=repAdjust[x],detectionThresh=detectionThresh))
			
			} else if(!is.null(meanD)&is.null(meanO)) {
		    	# single and double stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
		    	        .Call(.cpp.probabilityPeaksSD,
		    	            genotypeArray=cons$genotypes,
		                    alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                    heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                    DNAcont=rep(DNAcont,each=2), gradientS=locusGradient,
		                    meanD=meanD,
		                    meanS=locusStutter,
		                    scale=scale,degradation=rep(1+degradation,each=2),
		                    fragLengths=hypothesis$alleleDb[,2],
		                    fragNames=as.numeric(rownames(hypothesis$alleleDb)),
		                    stutterIndex = hypothesis$alleleDb[,3],
		                    repAdjust=repAdjust[x],detectionThresh=detectionThresh))

			} else if(is.null(meanD)&!is.null(meanO)) {
		    	# single and over stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
		    	        .Call(.cpp.probabilityPeaksSO,
		    	            genotypeArray=cons$genotypes,
		                    alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                    heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                    DNAcont=rep(DNAcont,each=2), gradientS=locusGradient,
		                    meanO=meanO,
		                    meanS=locusStutter,
		                    scale=scale,degradation=rep(1+degradation,each=2),
		                    fragLengths=hypothesis$alleleDb[,2],
		                    fragNames=as.numeric(rownames(hypothesis$alleleDb)),
			                stutterIndex = hypothesis$alleleDb[,3],
		                    repAdjust=repAdjust[x],detectionThresh=detectionThresh))
			}
		}
    # diagnose
	if(diagnose==TRUE) return(probs)
	# set probability of impossible genotypes to 0
	probs[is.na(probs)] = 0
        }
    return(probs)
    }



# get probability pf each genotype combination
peak.heights.per.locus = function(genotypeArray,alleles,
            heights,DNAcont,gradientS,meanD=NULL,
            meanO=NULL,meanS,scale,degradation,
            fragLengths,LUSvals,repAdjust=NULL,
            detectionThresh,diagnose=FALSE)
	{
	# get combined probability of all peaks for each genotype combination seperately
	Probs = apply(genotypeArray,MARGIN=2,FUN=function(x) probability.peaks(genotype=x,
	        alleles=alleles,heights=heights,
        	DNAcont=DNAcont,
        	gradientS=gradientS,meanD=meanD,
        	meanO=meanO,meanS=meanS,
        	scale=scale,degradation=degradation,
        	fragLengths=fragLengths,LUSvals=LUSvals,
        	repAdjust=repAdjust,
        	detectionThresh=detectionThresh,
	        diagnose=diagnose))
    # diagnose
	if(diagnose==TRUE) return(Probs)
	# give impossible values a probability of 0
	Probs = unlist(Probs)
	Probs[is.na(Probs)] = 0
	return(unlist(Probs))
	}


# get probability of every peak given current parameters
probability.peaks = function(genotype,alleles,
            heights,DNAcont,gradientS,meanD=NULL,
            meanO=NULL,meanS,scale,degradation,
            fragLengths,LUSvals,repAdjust=NULL,
            detectionThresh,diagnose=FALSE)
	{
	# convert genotype to numeric
	genotype = as.numeric(genotype)
	# get mean expected peak heights
	gammaMu = peak.height.dose(genotype=genotype,
	        alleles=alleles,heights=heights,
	        DNAcont=DNAcont,gradientS=gradientS,
	        meanD=meanD,meanO=meanO,meanS=meanS,
	        degradation=degradation,fragLengths=fragLengths,
	        LUSvals=LUSvals,repAdjust=repAdjust)
	names(heights) = alleles
	# give peak heights to dropout alleles (height=0)
	peakHeights = unlist(heights)
	gammaMus = gammaMu
	dropoutIndex = which(!names(gammaMus)%in%names(peakHeights))
	if(length(dropoutIndex)!=0)
	    {
        allelesToAdd = unique(names(gammaMus)[dropoutIndex])
        toAdd = rep(0,times=length(allelesToAdd))
        names(toAdd) = allelesToAdd
        peakHeights = c(peakHeights,toAdd)
        peakHeights = peakHeights[order(names(peakHeights))]
	    }
	# sort
	gammaMus = gammaMus[order(names(gammaMus))]
	peakHeights = peakHeights[order(names(peakHeights))]
	# scale = variance/mu
	# scale = (stanDev^2)/mu
	# alpha = (mu^2)/variance
	# alpha = (mu^2)/(stanDev^2)
	# shape = mu/scale
	# create shape and scale
	shapesVec = gammaMus/scale
	scalesVec = rep(scale,times=length(shapesVec)) 
    # diagnose
	if(diagnose==TRUE) return(list(height=peakHeights,
	mu=gammaMus,sigma=sqrt(shapesVec*scalesVec^2)))
	# probability densities
	pdf = vector(length=length(peakHeights))
	dropoutIndex = which(peakHeights==0)
	if(length(dropoutIndex)>0)
		{
		if(length(dropoutIndex)!=length(peakHeights))
		    {
		    # non-dropout = gamma point density
		    #pdf[-dropoutIndex] = mapply(FUN=function(x,k,t) dgamma(x=x,shape=k,scale=t), 
		    #        x=peakHeights[-dropoutIndex], 
		    #        k=shapesVec[-dropoutIndex], 
		    #        t=scalesVec[-dropoutIndex])
		    # discrete approximation to pmf
		    pdf[-dropoutIndex] = mapply(FUN=function(x,k,t) pgamma(q=x+0.5,shape=k,scale=t)-
		        pgamma(q=x-0.5,shape=k,scale=t), 
		            x=peakHeights[-dropoutIndex], 
		            k=shapesVec[-dropoutIndex], 
		            t=scalesVec[-dropoutIndex])
		            
		    }
		# dropout = integral from 0 to threshold
		pdf[dropoutIndex] = mapply(FUN=function(k,t) pgamma(q=detectionThresh,shape=k,scale=t), 
		        k=shapesVec[dropoutIndex], 
		        t=scalesVec[dropoutIndex])
		} else {
		# non-dropout = gamma point density
		pdf = mapply(FUN=function(x,k,t) pgamma(q=x+0.5,shape=k,scale=t)-
		    pgamma(q=x-0.5,shape=k,scale=t), 
		        x=peakHeights, 
		        k=shapesVec, 
		        t=scalesVec)
		}

	# set impossible values to 0 likelihood
	pdf[which(is.infinite(pdf))] = 0
	# output
	return(prod(pdf))
	}


# get mean expected peak height at each position
peak.height.dose = function(genotype,alleles,heights,
            DNAcont,gradientS,meanD=NULL,meanO=NULL,
            meanS,degradation,fragLengths,LUSvals,
            repAdjust=NULL)
	{
	# positions of stutter alleles
	stutterPos = genotype-1
	allPos = c(genotype,stutterPos)
	if(!is.null(meanD)) 
		{
		# positions of double stutter alleles
		doubleStutterPos = genotype-2
		allPos = c(allPos,doubleStutterPos)
		}
	if(!is.null(meanO))
		{
		# positions of over stutter alleles
		overStutterPos = genotype+1
		allPos = c(allPos,overStutterPos)
		}
	allPos = unique(round(allPos,1))
    # index frag lengths and LUS values
	fragLengthIndex = sapply(genotype,FUN=function(x) which(names(fragLengths)==x))
    # base dose
	dose = repAdjust*rep(DNAcont,each=2)*rep(1+degradation,each=2)^fragLengths[fragLengthIndex]
	# stutter rate
	stutterRate = meanS+(gradientS*LUSvals[fragLengthIndex])
    # dose from stutters
	muS = dose * stutterRate	
	if(!is.null(meanD))
		{
		# dose from double stutters
		muSd = dose * meanD
		stutterRate = stutterRate + meanD
		}
	if(!is.null(meanO))
		{
		# dose from over stutters
		muSo = dose * meanO
		stutterRate = stutterRate + meanO
		}
    # dose from allelic
	muA = dose * (1 - stutterRate)
	# tidy
	names(muA) = genotype
	names(muS) = stutterPos
	allPos = round(allPos,1)
	stutterPos = round(stutterPos,1)
	genotype = round(genotype,1)
	if(!is.null(meanD))
		{
		names(muSd) = doubleStutterPos
		doubleStutterPos = round(doubleStutterPos,1)
		}
	if(!is.null(meanO))
		{
		names(muSo) = overStutterPos
		overStutterPos = round(overStutterPos,1)
		}
    # combine doses for each position
	if(!is.null(meanD)&!is.null(meanO))
		{
		# S+D+O
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSd[which(doubleStutterPos==x)])+
		                                    sum(muSo[which(overStutterPos==x)]))
		names(muX) = allPos
		} else if(is.null(meanD)&!is.null(meanO)) {
		# S+O
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSo[which(overStutterPos==x)]))
		names(muX) = allPos
		} else if(!is.null(meanD)&is.null(meanO)) {
		# S+D
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)])+
		                                    sum(muSd[which(doubleStutterPos==x)]))
		names(muX) = allPos
		} else if(is.null(meanO)&is.null(meanD)) {
		# S
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+
		                                    sum(muS[which(stutterPos==x)]))
		names(muX) = allPos
		}
	return(muX)
	} 




# Penalties to apply to the likelihood.
# Documentation is in man directory.
penalties.peaks <- function(nloc, degradation=NULL,
                       degradationPenalty=50, gradientS = NULL, 
                       gradientAdjust=NULL,stutterAdjust=NULL,
                       stutterPenalty = 0.2,# stutterSD=0.2, 
                       meanS=NULL,meanD=NULL,meanO=NULL,
                       scale=NULL, scaleSD=1, ...) {
    # set initial penalty
    result = 1
    # Normalizes by number of loci so product of penalties same as in old code.
    # Some penalties are per locus and other for all locus. Hence this bit of run
    # around.
    normalization = 1.0 / max(nloc, 1)
    # penalty on degradation
    result = result * exp(-sum(degradation) * degradationPenalty 
                           * normalization )
    # penalty on gradientAdjust
    result = result * dnorm(log10(gradientAdjust),mean=0, sd=0.15)
    # penalty on stutterAdjust
    result = result * dnorm(log10(stutterAdjust),mean=0, sd=0.5)
    # penalty on meanD
    if(!missing(meanD) & !is.null(meanD))
        {
        # mean = 0.01, sd = 5e-5
        result = result * dgamma(meanD,shape=0.02/0.018,scale=0.018)
        }
    # penalty on meanO
    if(!missing(meanO) & !is.null(meanO))
	    {
	    # mean = 0.01, sd = 5e-5
	    result = result * dgamma(meanO,shape=0.02/0.018,scale=0.018)
	    }
  return(result)
}










    
