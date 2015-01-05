


# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.vectors.peaks <- function(hypothesis, addAttr=FALSE, likeMatrix=FALSE, diagnose=FALSE, ...) {

  hypothesis = likeLTD:::add.args.to.hypothesis(hypothesis, ...)

  sanity.check.peaks(hypothesis) # makes sure hypothesis has right type.
  locusCentric = transform.to.locus.centric.peaks(hypothesis)

  functions <- mapply(likeLTD:::create.likelihood.per.locus.peaks, locusCentric,
                      MoreArgs=list(addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose))

likeLTD:::create.likelihood.per.locus.peaks(locusCentric[[1]],addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose)


  #if(is.null(hypothesis$dropinPenalty)) hypothesis$dropinPenalty = 2 
  if(is.null(hypothesis$degradationPenalty)) hypothesis$degradationPenalty = 50 
  if(is.null(hypothesis$stutterPenalty)) hypothesis$stutterPenalty = 0.2


  likelihood.vectors <- function(degradation=NULL, DNAcont=NULL, scale=NULL, stutterMean=NULL, stutterAdjust=NULL, dropin=NULL, stutterGradient=NULL, 
                                doubleStutterRate = NULL,overStutterRate = NULL,
                                 repAdjust=NULL, detectionThresh=hypothesis$detectionThresh, degradationPenalty=hypothesis$degradationPenalty,stutterPenalty=hypothesis$stutterPenalty, ...) {
    # Call each and every function in the array.
    arguments = list(degradation=degradation, DNAcont=DNAcont,
                     scale = scale, repAdjust=repAdjust,
stutterMean=stutterMean,
stutterGradient=stutterGradient, 
doubleStutterRate = doubleStutterRate,overStutterRate=overStutterRate,
detectionThresh = detectionThresh,
                     degradationPenalty=degradationPenalty, stutterPenalty=stutterPenalty, dropin=dropin)
    callme <- function(objective,stut) {
    args = append(arguments, list(stutterAdjust=stut))
      do.call(objective, args)
    }
    if(length(stutterAdjust) == 1) stutterAdjust = rep(stutterAdjust, length(functions))
#    if(setequal(names(stutter), colnames(hypothesis$cspProfile)))
#      stutterAdjust <- stutterAdjust[colnames(hypothesis$cspProfile)]
    objectives = mapply(callme, functions, stutterAdjust)
    arguments = append(arguments, list(...))
if(diagnose==TRUE) return(objectives)
    # calculate penalties
    pens <- do.call(penalties.peaks, append(arguments,list(stutterAdjust=stutterAdjust,nloc=ncol(hypothesis$queriedProfile))))
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

  cons = likeLTD:::likelihood.constructs.per.locus.peaks(hypothesis)
  doR = !is.null(hypothesis$doR) && hypothesis$doR == TRUE

  result.function <- function(scale,stutterMean,stutterAdjust,
stutterGradient,
doubleStutterRate=NULL,overStutterRate=NULL,
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
    #   ...: Any other parameter, e.g. for penalty functions. 
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
	repRes <- likeLTD:::peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
				scale=scale, stutterMean=stutterMean, stutterAdjust=stutterAdjust, stutterGradient=stutterGradient,
				doubleStutterRate = doubleStutterRate, overStutterRate=overStutterRate,
degradation=degradation, 
				repAdjust=repAdjust,detectionThresh=detectionThresh,doR=doR,diagnose=diagnose)
	return(repRes)
	}

    repRes <- matrix(likeLTD:::peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
				scale=scale, stutterMean=stutterMean,stutterAdjust=stutterAdjust, stutterGradient=stutterGradient, 
				doubleStutterRate = doubleStutterRate,overStutterRate=overStutterRate,
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
  HYPOTHESIS <<- hypothesis

genotypes = likeLTD:::explain.all.peaks(cspPresence,knownPresence,hypothesis$knownProfs,alleles,hypothesis$nUnknowns,hypothesis$peaksProfile,hypothesis$heightsProfile,hypothesis$doDoubleStutter,hypothesis$doOverStutter)

# get index of which alleles are from known contributors - do not want population allele probabilities for known contributors
if(nrow(hypothesis$knownProfs)>0) 
	{
	kIndex = (hypothesis$nUnknowns*2)+(1:(2*nrow(hypothesis$knownProfs)))
	} else {
	kIndex = c()
	}

GENS <<- genotypes
KINDEX <<- kIndex

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
peaks.probabilities = function(hypothesis,cons,DNAcont,scale,stutterMean,stutterAdjust,stutterGradient,
       doubleStutterRate=NULL,overStutterRate=NULL,
degradation,repAdjust,detectionThresh,doR=FALSE,diagnose=FALSE)#,doC=TRUE)
    {
HYPOTHESIS<<-hypothesis
CONS<<-cons
DNACONT<<-DNAcont
SCALE<<-scale
STUTTERMEAN<<-stutterMean
STUTTERADJUST<<-stutterAdjust
STUTTERGRADIENT<<-stutterGradient
DOUBLESTUTTERRATE<<-doubleStutterRate
OVERSTUTTERRATE<<-overStutterRate
DEGRADATION<<-degradation
REPADJUST<<-repAdjust
DETECTIONTHRESH<<-detectionThresh
    # return a function that computes the 
    if(hypothesis$doDropin==TRUE)
        {
        stop("Cannot handle dropin")
        } else {
    #genotypes = matrix(as.numeric(rownames(hypothesis$alleleDb))[cons$genotypes],ncol=ncol(cons$genotypes))
	#genotypes = matrix(as.numeric(genotypes),ncol=ncol(genotypes))

	
CONS <<- cons
#	GENOTYPES <<- genotypes
#	PEAKS <<- as.numeric(hypothesis$peaksProfile[[1]])
#	HEIGHTS <<- unlist(as.numeric(hypothesis$heightsProfile[[1]]))
#	SIZES <<- unlist(hypothesis$sizesProfile[[1]])
#	DNACONT <<- DNAcont
#	STUTTERMEAN <<- stutterMean
#	STUTTERADJUST <<- stutterAdjust
#	SCALE <<- scale
#	DEG <<- degradation
#	FRAGLENGTHS <<- hypothesis$alleleDb[,2]
#	FRAGNAMES <<- as.numeric(rownames(hypothesis$alleleDb))
#	REPADJUST <<- repAdjust[1]
#	THRESHOLD <<- threshold
DSR<<- doubleStutterRate
	if(doR==TRUE|diagnose==TRUE)
		{
        	# probabilities for each replicate
        	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) peak.heights.per.locus(genotypeArray=cons$genotypes,
										alleles=hypothesis$peaksProfile[[x]],heights=hypothesis$heightsProfile[[x]],
										sizes=hypothesis$sizesProfile[[x]],DNAcont=DNAcont,
										stutterMean=stutterMean,stutterAdjust=stutterAdjust,
										doubleStutterRate=doubleStutterRate,
										overStutterRate=overStutterRate,
										scale=scale,degradation=degradation,
										fragLengths=hypothesis$alleleDb[,2],repAdjust=repAdjust[x],
										detectionThresh=detectionThresh,diagnose=diagnose))
		} else {
		if(!is.null(doubleStutterRate)&!is.null(overStutterRate))
			{
			# single, double and over stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) .Call(.cpp.probabilityPeaksSDO,genotypeArray=cons$genotypes,
		                                                                        alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                                                                        heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                                                                        DNAcont=rep(DNAcont,each=2),stutterMean=stutterMean, stutterAdjust=stutterAdjust,
		                                                                        stutterGradient=stutterGradient,
		                                                                        doubleStutterRate=doubleStutterRate,overStutterRate=overStutterRate,
		                                                                        scale=scale,degradation=rep(1+degradation,each=2),
		                                                                        fragLengths=hypothesis$alleleDb[,2],
		                                                                        fragNames=as.numeric(rownames(hypothesis$alleleDb)),
		                                                                        repAdjust=repAdjust[x],detectionThresh=detectionThresh))

			} else if(is.null(doubleStutterRate)&is.null(overStutterRate)) {
		   	# single stutter only
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) .Call(.cpp.probabilityPeaksS,genotypeArray=cons$genotypes,
		                                                                        alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                                                                        heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                                                                        DNAcont=rep(DNAcont,each=2),stutterMean=stutterMean,stutterAdjust=stutterAdjust,
		                                                                        stutterGradient=stutterGradient,
		                                                                        scale=scale,degradation=rep(1+degradation,each=2),
		                                                                        fragLengths=hypothesis$alleleDb[,2],
		                                                                        fragNames=as.numeric(rownames(hypothesis$alleleDb)),
		                                                                        repAdjust=repAdjust[x],detectionThresh=detectionThresh))
			
			} else if(!is.null(doubleStutterRate)&is.null(overStutterRate)) {
		    	# single and double stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) .Call(.cpp.probabilityPeaksSD,genotypeArray=cons$genotypes,
		                                                                        alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                                                                        heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                                                                        DNAcont=rep(DNAcont,each=2),stutterMean=stutterMean, stutterAdjust=stutterAdjust,
		                                                                        stutterGradient=stutterGradient,
		                                                                        doubleStutterRate=doubleStutterRate,
		                                                                        scale=scale,degradation=rep(1+degradation,each=2),
		                                                                        fragLengths=hypothesis$alleleDb[,2],
		                                                                        fragNames=as.numeric(rownames(hypothesis$alleleDb)),
		                                                                        repAdjust=repAdjust[x],detectionThresh=detectionThresh))

			} else if(is.null(doubleStutterRate)&!is.null(overStutterRate)) {
		    	# single and over stutter
		    	probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) .Call(.cpp.probabilityPeaksSO,genotypeArray=cons$genotypes,
		                                                                        alleles=as.numeric(hypothesis$peaksProfile[[x]]),
		                                                                        heights=unlist(as.numeric(hypothesis$heightsProfile[[x]])),
		                                                                        DNAcont=rep(DNAcont,each=2),stutterMean=stutterMean, stutterAdjust=stutterAdjust,
		                                                                        stutterGradient=stutterGradient,
		                                                                        overStutterRate=overStutterRate,
		                                                                        scale=scale,degradation=rep(1+degradation,each=2),
		                                                                        fragLengths=hypothesis$alleleDb[,2],
		                                                                        fragNames=as.numeric(rownames(hypothesis$alleleDb)),
		                                                                        repAdjust=repAdjust[x],detectionThresh=detectionThresh))
			}
		}
	if(diagnose==TRUE) return(probs)
	probs[is.na(probs)] = 0
        }
    return(probs)
    }



# Function to run peak heights per locus
# Parameters:
# genotypeArray = Matrix of genotype combinations for this locus
# alleles = CSP alleles for this locus
# heights = CSP peak heights for this locus
# sizes = CSP allele sizes for this locus
# rcont = current rcont value
# rhoA = allelic constant parameter, single value
# rhoS = stutter constant parameter, single value
peak.heights.per.locus = function(genotypeArray,alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,doubleStutterRate=NULL,overStutterRate=NULL,scale,degradation,fragLengths,repAdjust=NULL,detectionThresh,diagnose=FALSE)
	{
	#index = !is.na(alleles)
	#alleles = alleles[index]
	#heights = heights[index]
	#sizes = sizes[index]


	# result vector
	#if(parallel==FALSE)
	#	{
	        Probs = apply(genotypeArray,MARGIN=2,FUN=function(x) probability.peaks(x,alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,doubleStutterRate,scale,degradation,fragLengths,repAdjust,detectionThresh,diagnose))
	#	} else {
	#        Probs = mclapply(1:ncol(genotypeArray),FUN=function(x) probability.peaks(genotypeArray[,x],alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,scale,degradation,fragLengths,repAdjust,diagnose),mc.cores=cores)
	#	}
	if(diagnose==TRUE) return(Probs)
	Probs = unlist(Probs)
	Probs[is.na(Probs)] = 0
	return(unlist(Probs))
	}


# Function to get probability of observed peak heights for a 
# single genotype combination when there is no dropout required
# to explain the genotype combination
# Parameters:
# genotype = single genotype combination (vector of 2*NU)
# alleles = CSP alleles for this locus
# heights = CSP peak heights for this locus
# sizes = CSP allele sizes for this locus
# rcont = current rcont value
# rhoA = allelic constant parameter, single value
# rhoS = stutter constant parameter, single value
probability.peaks = function(genotype,alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,doubleStutterRate=NULL,overStutterRate=NULL,scale,degradation,fragLengths,repAdjust=NULL,detectionThresh,diagnose=FALSE)
	{	
	genotype = as.numeric(genotype)

	# get means
	gammaMu = peak.height.dose(genotype,alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,doubleStutterRate,overStutterRate,degradation,fragLengths,repAdjust)
	names(heights) = alleles
	# give peak heights to dropout alleles
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
	# give means to peaks that are not hypothesised (all means should be zero)
	unhypIndex = which(!names(peakHeights)%in%names(gammaMus))
	if(length(unhypIndex)!=0)
		{
		toAdd = rep(0,times=length(unhypIndex))
		names(toAdd) = names(peakHeights)[unhypIndex]
		gammaMus = c(gammaMus, toAdd)
		}
	# sort
	gammaMus = gammaMus[order(names(gammaMus))]
	peakHeights = peakHeights[order(names(peakHeights))]

	# scale = variance/mu
	# scale = (stanDev^2)/mu
	# alpha = (mu^2)/variance
	# alpha = (mu^2)/(stanDev^2)

	# shape = mu/scale
	shapesVec = gammaMus/scale
	scalesVec = rep(scale,times=length(shapesVec)) 

	

	if(diagnose==TRUE) return(list(height=peakHeights,mu=gammaMus,sigma=sqrt(shapesVec*scalesVec^2)))#shapes=shapesVec,scales=scalesVec))

	# probability densities
	pdf = vector(length=length(peakHeights))
	dropoutIndex = which(peakHeights==0)

	if(length(dropoutIndex)>0)
		{
		# non-dropout = gamma point density
		pdf[-dropoutIndex] = mapply(FUN=function(x,k,t) dgamma(x=x,shape=k,scale=t), x=peakHeights[-dropoutIndex], k=shapesVec[-dropoutIndex], t=scalesVec[-dropoutIndex])
		# dropout = integral from 0 to threshold
		pdf[dropoutIndex] = mapply(FUN=function(k,t) pgamma(q=detectionThresh,shape=k,scale=t), k=shapesVec[dropoutIndex], t=scalesVec[dropoutIndex])
		} else {
		pdf = mapply(FUN=function(x,k,t) dgamma(x=x,shape=k,scale=t), x=peakHeights, k=shapesVec, t=scalesVec)
		}



	# set impossible values to 0 likelihood
	pdf[which(is.infinite(pdf))] = 0
	# output
	return(prod(pdf))
	}


# get "dose" (alpha values) for a single locus and single genotype
# dropout alleles are given an alpha value of 0
# Parameter:
# genotype = genotypes of hypothesised contributors  (from genotype matrix)
# alleles = alleles observed in CSP (including stutters etc)
# heights = peak heights from CSP
# sizes = base pair sizes of CSP peaks
# DNAcont = relative DNA contribution of contributors
# DNAproxy = Assumed total amount of DNA in the sample (usually the sum of peak heights at this locus)
# rhoA = allelic constant parameter
# rhoS = stutter constant parameter
# scale = stanDev constant parameter for gamma
peak.height.dose = function(genotype,alleles,heights,sizes,DNAcont,stutterMean,stutterAdjust,doubleStutterRate=NULL,overStutterRate=NULL,degradation,fragLengths,repAdjust=NULL)
	{
	# positions of stutter alleles
	stutterPos = genotype-1
	allPos = c(genotype,stutterPos)
	if(!is.null(doubleStutterRate)) 
		{
		doubleStutterPos = genotype-2
		allPos = c(allPos,doubleStutterPos)
		}
	if(!null(overStutterRate))
		{
		overStutterPos = genotype+1
		allPos = c(allPos,overStutterPos)
		}
	allPos = unique(round(allPos,1))
	# all positions
	# round is needed because of floating point error
	# e.g. 31.2 != 32.2-1
	#allPos = sort(unique(round(c(stutterPos,genotype),1)))
	# add sizes of genotyp alleles that are not in the CSP (dropout alleles)
	names(sizes) = alleles
	#allSizes = sizes
	#if(any(!genotype%in%alleles)) allSizes = fillSizes(sizes,alleles,genotype,allPos)
	# changes sizes to order of genotype rather than order of CSP
	#sizesGen = sapply(genotype,FUN=function(x) allSizes[which(names(allSizes)==x)])
	fragLengthIndex = sapply(genotype,FUN=function(x) which(names(fragLengths)==x))
    	# degradation adjustment
	degAdjust = repAdjust*rep(DNAcont,each=2)*rep(1+degradation,each=2)^fragLengths[fragLengthIndex]
#return(DNAcont)
	#degAdjust = degAdjust/sum(degAdjust)   # this was previously omega - converts rcont to proportion
	stutterRate = stutterMean*stutterAdjust
	# stutter alpha
	#alphaS = degAdjust * rhoS * DNAproxy / sizesGen
	muS = degAdjust * stutterRate #* DNAproxy / sizesGen
	if(!is.null(doubleStutterRate))
		{
		muSd = degAdjust * doubleStutterRate
		stutterRate = stutterRate + doubleStutterRate
		}
	if(!is.null(overStutterRate))
		{
		muSo = degAdjust * overStutterRate
		stutterRate = stutterRate + overStutterRate
		}
	# allelic alpha
	#alphaA = degAdjust * rhoA * DNAproxy / sizesGen
	muA = degAdjust * (1 - stutterRate) #* DNAproxy / sizesGen
	names(muA) = genotype
	names(muS) = stutterPos
	allPos = round(allPos,1)
	stutterPos = round(stutterPos,1)
	genotype = round(genotype,1)
	if(!is.null(doubleStutterRate))
		{
		names(muSd) = doubleStutterPos
		doubleStutterPos = round(doubleStutterPos,1)
		}
	if(!is.null(overStutterRate))
		{
		names(muSo) = overStutterPos
		overStutterPos = round(overStutterPos,1)
		}

#return(fragLengthIndex-1)
	# get total alphas for each position
	if(!is.null(doubleStutterRate)&!is.null(overStutterRate))
		{
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+sum(muS[which(stutterPos==x)])+sum(muSd[which(doubleStutterPos==x)])+sum(muSo[which(overStutterPos==x)]))
		} else if(is.null(doubleStutterRate)&!is.null(overStutterRate)) {
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+sum(muS[which(stutterPos==x)])+sum(muSo[which(overStutterPos==x)]))
		names(muX) = allPos
		} else if(!is.null(doubleStutterRate)&is.null(overStutterRate)) {
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+sum(muS[which(stutterPos==x)])+sum(muSd[which(doubleStutterPos==x)]))
		names(muX) = allPos
		} else if(is.null(overStutterRate)&is.null(doubleStutterRate)) {
		muX = sapply(allPos,FUN=function(x) sum(muA[which(genotype==x)])+sum(muS[which(stutterPos==x)]))
		names(muX) = allPos
		}
	return(muX)
	} 


# function to add sizes for positions in the hypothesised
# genotype that were not observed in the CSP
# useful when there is stutter and/or dropout
# Parameter:
# sizes = base pair sizes of CSP peaks
# alleles = alleles observed in CSP (including stutters etc)
# genotype = genotypes of hypothesised contributors (from genotype matrix)
# allPos = allelic positions including stutter and dropout etc
fillSizes = function(sizes,alleles,genotype,allPos)
	{
	# which alleles are missing from sizes
	index = unique(allPos[which(!allPos%in%names(sizes))])
	newsizes = unlist(sizes)
	for(i in 1:length(index)) newsizes = addMissingAlleleSize(index[i],newsizes)
	return(newsizes[order(names(newsizes))])
	}

# function to add sizes for missing alleles based on their allele name difference
# from the nearest allele with a size
# e.g. the size of allele 16 is the size of allele 18 - (2*4)
# need to fix this for .2 if nearet is not .2
addMissingAlleleSize = function(index,sizes)
	{
	# size difference between missing alleles and sizes
	diffSizes = index - as.numeric(names(sizes))
	# get which size allele is closest to the missing allele
	closestIndex = which.min(abs(diffSizes))
	# get the new size
	newSizes = sizes[closestIndex] + (4*diffSizes[closestIndex]) 
	# add to sizes
	outsizes = c(newSizes,sizes)
	names(outsizes)[1] = index
	return(outsizes)
	}









# Penalties to apply to the likelihood.
# Documentation is in man directory.
penalties.peaks <- function(nloc, degradation=NULL,
                       degradationPenalty=50, stutterAdjust=NULL,
                       stutterPenalty = 0.2,# stutterSD=0.2, 
stutterMean=NULL,doubleStutterRate=NULL,overStutterRate=NULL,
scale=NULL, scaleSD=1, ...) {
  result = 1
  # Normalizes by number of loci so product of penalties same as in old code.
  # Some penalties are per locus and other for all locus. Hence this bit of run
  # around.
  normalization = 1.0 / max(nloc, 1)

  if(!missing(degradation) & !is.null(degradation))
    result = result * exp(-sum(degradation) * degradationPenalty 
                           * normalization )

    # gaussian penalty on stutter adjustment
    # sd of stutter percentage ranges from 1.2 to 3.2 (will set penalty based on 3 here)
    #(see Leclair-et-al (2004) Systematic Analysis of Stutter Percentages and Allele Peak Height and Peak Area Ratios at Heterozygous STR Loci for Forensic Casework and Database Samples)
    if(!missing(stutterAdjust) & !is.null(stutterAdjust))
        {
        result = result * dnorm(log10(stutterAdjust),mean=0, sd=stutterPenalty)
        }

 if(!missing(doubleStutterRate) & !is.null(doubleStutterRate))
	{
	#result = result * 10^(log10(alteration)/10)
	result = result * dnorm(log(doubleStutterRate),mean=-14, sd=4.6)
	}

 if(!missing(overStutterRate) & !is.null(overStutterRate))
	{
	#result = result * 10^(log10(alteration)/10)
	result = result * dnorm(log(overStutterRate),mean=-14, sd=4.6)
	}


  return(result)
}










    
