# Creates a likelihood function from the input hypothesis
# Documentation is in man directory.
create.likelihood.vectors.peaks <- function(hypothesis, addAttr=FALSE, 
				likeMatrix=FALSE, diagnose=FALSE, ...) 
	{
	# add arguments to hypothesis
	hypothesis = add.args.to.hypothesis(hypothesis, ...)
	# check hypothesis has the right type
	sanity.check.peaks(hypothesis)
	# convert hypothesis to locus specific
	locusCentric = transform.to.locus.centric.peaks(hypothesis)
	# functions to perform on each locus
	functions <- mapply(create.likelihood.per.locus.peaks, locusCentric,
                      MoreArgs=list(addAttr=addAttr, likeMatrix=likeMatrix, diagnose=diagnose))
	#create.likelihood.per.locus.peaks(locusCentric[[1]],addAttr=addAttr, 
	#				likeMatrix=likeMatrix, diagnose=diagnose)
	# set degradation penalty if it is missing
	if(is.null(hypothesis$degradationPenalty)) hypothesis$degradationPenalty = 50 
	# output function, to be run every iteration
	likelihood.vectors <- function(degradation=NULL, DNAcont=NULL, 
					scale=NULL, gradientS = NULL, 
					gradientAdjust=NULL, interceptAdjust=NULL, 
					dropin=NULL, interceptS=NULL, 
					meanD=NULL, meanO=NULL, repAdjust=NULL, 
					detectionThresh=hypothesis$detectionThresh, 
					degradationPenalty=hypothesis$degradationPenalty,
					stutterPenalty=hypothesis$stutterPenalty, ...) 
		{
    		# set arguments
		arguments = list(degradation=degradation, DNAcont=DNAcont,
				scale = scale, repAdjust=repAdjust,
				gradientS = gradientS,interceptS=interceptS, 
				meanD = meanD,meanO=meanO,
				detectionThresh = detectionThresh,
		                degradationPenalty=degradationPenalty, 
				stutterPenalty=stutterPenalty, dropin=dropin)
		# single locus objective function
		callme <- function(objective,grad,stut) 
			{
			args = append(arguments, list(gradientAdjust=grad,interceptAdjust=stut))
			do.call(objective, args)
			}
		if(length(gradientAdjust) == 1) gradientAdjust = rep(gradientAdjust, length(functions))
		if(length(interceptAdjust) == 1) interceptAdjust = rep(interceptAdjust, length(functions))
		# call every objective function (one per locus)
		objectives = mapply(callme, functions, gradientAdjust,interceptAdjust)
		arguments = append(arguments, list(...))
		if(diagnose==TRUE) return(objectives)
		# calculate penalties
		pens <- do.call(penalties.peaks, append(arguments,list(gradientAdjust=gradientAdjust,
				interceptAdjust=interceptAdjust,nloc=ncol(hypothesis$queriedProfile))))
    		list(objectives=objectives, penalties=pens)
  		}
	if(addAttr) 
		{
		attr(likelihood.vectors, "hypothesis") <- hypothesis
		attr(likelihood.vectors, "functions") <- functions
		}
	return(likelihood.vectors)
	}


create.likelihood.per.locus.peaks <- function(hypothesis, addAttr=FALSE, 
				likeMatrix = FALSE, diagnose=FALSE) 
	{
	# Creates a likelihood function for a given hypothesis and locus
	#
	# A hypothesis is given by the number of unknown contributors, whether to model
	# dropin, so on and so forth.
	cons = likelihood.constructs.per.locus.peaks(hypothesis)
	doR = !is.null(hypothesis$doR) && hypothesis$doR == TRUE

	result.function <- function(scale,gradientS,gradientAdjust,interceptAdjust,
				interceptS,meanD=NULL,meanO=NULL,repAdjust=NULL,
				degradation=NULL, DNAcont=NULL, 
				detectionThresh = NULL, dropin=NULL, ...) 
		{
		# Likelihood function for a given hypothesis and locus
		#
		# This function is specific to the hypothesis for which it was created.
		# It hides everything except the nuisance parameters over which to
		# optimize.
		#
		# Parameters:
		#   degradation: relative degradation from each profiled individual in this
		#                hypothesis
	    	#   ...: Any other parameter, e.g. for penalty functions. genotypeArray=cons$genotypes,
	    	#        These parameters are ignored here.
	    	# Returns: A scalar value giving the likelihood for this locus and
	    	#          hypothesis
		# complete repAdjust vector
		repAdjust = c(1,repAdjust)
		# perform some checks
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
			# diagnose result (for computing Z values)
			repRes <- peaks.probabilities(hypothesis=hypothesis, cons=cons, 
							DNAcont=DNAcont,scale=scale,
							gradientS = gradientS, 
							gradientAdjust=gradientAdjust,
							interceptAdjust=interceptAdjust, 
							interceptS=interceptS,
							meanD = meanD, meanO=meanO,
							degradation=degradation, 
							repAdjust=repAdjust,
							detectionThresh=detectionThresh,
							dropin=dropin,
							doR=doR,diagnose=diagnose)
			return(repRes)
			}
		# result
		res <- peaks.probabilities(hypothesis=hypothesis, cons=cons, DNAcont=DNAcont, 
					scale=scale,gradientS = gradientS,gradientAdjust=gradientAdjust, 						interceptAdjust=interceptAdjust,interceptS=interceptS, 
					meanD = meanD,meanO=meanO,degradation=degradation, 
					repAdjust=repAdjust,dropin=dropin,
					detectionThresh=detectionThresh,doR=doR)
		# multiply by allele probabilities
		factorsRes = res*cons$factors
		# Figure out likelihood for good and return.
		if(likeMatrix==FALSE)
			{
			return(sum(factorsRes))
			} else {
			return(factorsRes)
			}
		}
	# add some attributes
	if(addAttr) 
		{
		attr(result.function, "hypothesis") <- hypothesis
		attr(result.function, "constructs") <- cons
		}
	# return the result function to be performed every iteration
	result.function
	}


get.database.values = function(alleleDb, doOverStutter, doDoubleStutter)
	{
	# convert database alleles to numeric
	db = as.numeric(rownames(alleleDb))
	# add stutter allele to db vals
	out = c(db, db-1)
	# add over stutter values
	if(doOverStutter) out = c(out, db+1)
	# add double stutter values
	if(doDoubleStutter) out = c(out, db-2)
	sort(unique(round(out,1)))
	}


likelihood.constructs.per.locus.peaks = function(hypothesis) 
	{
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
	if(!is.matrix(cspPresence)) cspPresence = matrix(ncol=0, nrow=length(alleles))
	if(!is.matrix(knownPresence)) knownPresence = matrix(ncol=0, nrow=length(alleles))
	missingReps = apply(hypothesis$binaryProfile, 1, is.na)
	# get intial genotype array
	genotypes = explain.all.peaks(cspPresence,knownPresence,hypothesis$knownProfs,
					alleles,hypothesis$nUnknowns,hypothesis$peaksProfile,
					hypothesis$heightsProfile,hypothesis$doDoubleStutter,
					hypothesis$doOverStutter,hypothesis$doDropin)
	# get database values, including all stutter positions
	dbVals = get.database.values(hypothesis$alleleDb,hypothesis$doOverStutter,
					hypothesis$doDoubleStutter)
	# get index of which alleles are from known contributors
	#do not want population allele probabilities for known contributors
	if(nrow(hypothesis$knownProfs)>0) 
		{
		kIndex = (hypothesis$nUnknowns*2)+(1:(2*nrow(hypothesis$knownProfs)))
		} else {
		kIndex = c()
		}
	# get genotype probabilities, taking into account relatedness
	if(length(kIndex>0)) 
		{
		# Only take into account relatedness between Q and X under Hd 
		if(hypothesis$hypothesis=="defence")
			{
			# exclude known genotypes
		  	factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], 
					hypothesis$alleleDb,hypothesis$nUnknowns, 
					hypothesis$doDropin,hypothesis$queriedProfile,
                             		hypothesis$relatedness) 
			} else {
			# exclude known genotypes
  			factors = genotype.factors(genotypes[-kIndex,,drop=FALSE], 
					hypothesis$alleleDb,hypothesis$nUnknowns, 
					hypothesis$doDropin,hypothesis$queriedProfile,
                             		c(0,0))
			}
		} else {
		# Only take into account relatedness between Q and X under Hd 
		if(hypothesis$hypothesis=="defence")
			{
			factors = genotype.factors(genotypes, hypothesis$alleleDb,
                             		hypothesis$nUnknowns, hypothesis$doDropin,
					hypothesis$queriedProfile,hypothesis$relatedness) 
			} else {
			factors = genotype.factors(genotypes, hypothesis$alleleDb,
					hypothesis$nUnknowns, hypothesis$doDropin,
					hypothesis$queriedProfile,c(0,0))
			}
		}
	# convert genotypes from indices to repeat number
	genotypes = matrix(as.numeric(rownames(hypothesis$alleleDb))[genotypes],ncol=ncol(genotypes))
	# output list
	list(cspPresence=cspPresence, knownPresence=knownPresence,
        	missingReps=missingReps,genotypes=genotypes, factors=factors,
       		freqMat=hypothesis$alleleDb[, 1],dbVals=dbVals)
	}


# function to be called at each iteration of maximisation
peaks.probabilities = function(hypothesis,cons,DNAcont,scale,
			gradientS,gradientAdjust,interceptAdjust,
			interceptS,meanD=NULL,meanO=NULL,degradation,
			repAdjust,detectionThresh,dropin=NULL,doR=FALSE,diagnose=FALSE)
	{
	# combine mean and adjustment
	locusGradient = gradientS*gradientAdjust
	locusIntercept = interceptS*interceptAdjust


TgenotypeArray<<-cons$genotypes
					TDNAcont<<-rep(DNAcont,each=2) 
					TgradientS<<-locusGradient
					TmeanD<<-meanD;TmeanO<<-meanO
					TinterceptS<<-locusIntercept
					Tdegradation<<-rep(1+degradation,each=2)
					TfragLengths<<-hypothesis$alleleDb[,2]
					TfragNames<<-as.numeric(rownames(hypothesis$alleleDb))
					TLUSvals <<- hypothesis$alleleDb[,3]
					Talleles<<-hypothesis$peaksProfile
					Theights<<-hypothesis$heightsProfile
					TrepAdjust<<-repAdjust;Tscale<<-scale
					TdetectionThresh<<-detectionThresh
					TdatabaseVals <<- cons$dbVals
					TfragProbs<<-hypothesis$alleleDb[,1]; Tdropin<<-dropin
	if(hypothesis$doDropin==TRUE)
		{
		# no dropin currently
		if(doR==TRUE|diagnose==TRUE)
			{
			# probabilities for each replicate
			probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
				peak.heights.per.locus(genotypeArray=cons$genotypes,
						alleles=hypothesis$peaksProfile[[x]],
						heights=hypothesis$heightsProfile[[x]],
						DNAcont=DNAcont,
						gradientS = locusGradient,
						meanD=meanD,meanO=meanO,
						interceptS=locusIntercept,
						scale=scale,degradation=degradation,
						fragLengths=hypothesis$alleleDb[,2],
						fragProbs=hypothesis$alleleDb[,1],
						LUSvals=hypothesis$alleleDb[,3],
						repAdjust=repAdjust[x],
						detectionThresh=detectionThresh,
						dropin=dropin,
						diagnose=diagnose))
			} else {
				# single, double and over stutter
		    		probs = .Call(.cpp.getProbabilitiesSDO_dropin,
					genotypeArray=cons$genotypes,
					DNAcont=rep(DNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals,
					fragProbs=hypothesis$alleleDb[,1], dropin=dropin)
			}
		} else {
		if(doR==TRUE|diagnose==TRUE)
			{
			# probabilities for each replicate
			probs = sapply(1:length(hypothesis$peaksProfile), FUN=function(x) 
				peak.heights.per.locus(genotypeArray=cons$genotypes,
						alleles=hypothesis$peaksProfile[[x]],
						heights=hypothesis$heightsProfile[[x]],
						DNAcont=DNAcont,
						gradientS = locusGradient,
						meanD=meanD,meanO=meanO,
						interceptS=locusIntercept,
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
		    		probs = .Call(.cpp.getProbabilitiesSDO,
					genotypeArray=cons$genotypes,
					DNAcont=rep(DNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
			} else if(is.null(meanD)&is.null(meanO)) {
		   		# single stutter only
		    		probs = .Call(.cpp.getProbabilitiesS,
					genotypeArray=cons$genotypes,
					DNAcont=rep(DNAcont,each=2), 
					gradientS=locusGradient,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
			
			} else if(!is.null(meanD)&is.null(meanO)) {
		    		# single and double stutter
		    		probs = .Call(.cpp.getProbabilitiesSD,
					genotypeArray=cons$genotypes,
					DNAcont=rep(DNAcont,each=2), 
					gradientS=locusGradient,
					meanD=meanD,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)

			} else if(is.null(meanD)&!is.null(meanO)) {
		    		# single and over stutter
		    		probs = .Call(.cpp.getProbabilitiesSO,
					genotypeArray=cons$genotypes,
					DNAcont=rep(DNAcont,each=2), 
					gradientS=locusGradient,
					meanO=meanO,
					interceptS=locusIntercept,
					degradation=rep(1+degradation,each=2),
					fragLengths=hypothesis$alleleDb[,2],
					fragNames=as.numeric(rownames(hypothesis$alleleDb)),
					LUSvals = hypothesis$alleleDb[,3],
					alleles=hypothesis$peaksProfile,
					heights=hypothesis$heightsProfile,
					repAdjust=repAdjust,scale=scale,
					detectionThresh=detectionThresh,
					databaseVals = cons$dbVals)
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
peak.heights.per.locus = function(genotypeArray,alleles,heights,DNAcont,
			gradientS,meanD=NULL,meanO=NULL,interceptS,scale,
			degradation,fragLengths,fragProbs=NULL,LUSvals,repAdjust=NULL,
			detectionThresh,dropin=FALSE,diagnose=FALSE)
	{
	# get combined probability of all peaks for each genotype combination seperately
	Probs = apply(genotypeArray,MARGIN=2,FUN=function(x) probability.peaks(genotype=x,
	        alleles=alleles,heights=heights,
        	DNAcont=DNAcont,
        	gradientS=gradientS,meanD=meanD,
        	meanO=meanO,interceptS=interceptS,
        	scale=scale,degradation=degradation,
        	fragLengths=fragLengths,fragProbs=fragProbs,LUSvals=LUSvals,
        	repAdjust=repAdjust,
        	detectionThresh=detectionThresh,dropin=dropin,
	        diagnose=diagnose))
	# diagnose
	if(diagnose==TRUE) return(Probs)
	# give impossible values a probability of 0
	Probs = unlist(Probs)
	Probs[is.na(Probs)] = 0
	return(unlist(Probs))
	}


# get probability of every peak given current parameters
probability.peaks = function(genotype,alleles,heights,DNAcont,
		gradientS,meanD=NULL,meanO=NULL,interceptS,scale,
		degradation,fragLengths,fragProbs=NULL,LUSvals,repAdjust=NULL,
		detectionThresh,dropin=NULL,diagnose=FALSE)
	{
	# convert genotype to numeric
	genotype = as.numeric(genotype)
	# get mean expected peak heights
	gammaMu = peak.height.dose(genotype=genotype,
	        alleles=alleles,heights=heights,
	        DNAcont=DNAcont,gradientS=gradientS,
	        meanD=meanD,meanO=meanO,interceptS=interceptS,
	        degradation=degradation,fragLengths=fragLengths,
	        LUSvals=LUSvals,repAdjust=repAdjust,dropin=dropin)
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
peak.height.dose = function(genotype,alleles,heights,DNAcont,
			gradientS,meanD=NULL,meanO=NULL,
            		interceptS,degradation,fragLengths,fragProbs=NULL,
			LUSvals,repAdjust=NULL,dropin=NULL)
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
	dose = repAdjust*rep(DNAcont,each=2)*rep(1+degradation,each=2)^-fragLengths[fragLengthIndex]
	# stutter rate
	stutterRate = interceptS+(gradientS*LUSvals[fragLengthIndex])
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
	# add dropin expected dose
	if(!is.null(fragProbs))
		{
		for(i in 1:length(muX))
			{
			index = which(round(as.numeric(names(fragProbs)),1)==names(muX)[i])
			if(length(index)>0)
				{
				muX[i] = muX[i]+fragProbs[index]*dropin 
				}
			}
		}
	return(muX)
	} 




# Penalties to apply to the likelihood.
# Documentation is in man directory.
penalties.peaks <- function(nloc, degradation=NULL,
                       degradationPenalty=50, gradientS = NULL, 
                       gradientAdjust=NULL,interceptAdjust=NULL,
                       stutterPenalty = 0.2,# stutterSD=0.2, 
                       interceptS=NULL,meanD=NULL,meanO=NULL,
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
    # penalty on gradientS
    # mean = 0.013, var=0.1 (actual var=0.01^2) 
    result = result * (dgamma(meanD,shape=0.013/(0.01^2/0.013),scale=0.01^2/0.013) * normalization)
    # penalty on interceptS
    # mean = 0.001, var=0.1 (actual mean=-0.076, actual var=0.13^2) 
    result = result * (dgamma(meanD,shape=0.001/(0.0001/0.001),scale=0.0001/0.001) * normalization)
    # penalty on gradientAdjust
    result = result * dnorm(log10(gradientAdjust),mean=0, sd=0.15)
    # penalty on interceptAdjust
    result = result * dnorm(log10(interceptAdjust),mean=0, sd=0.5)
    # penalty on meanD
    if(!missing(meanD) & !is.null(meanD))
        {
        # mean = 0.01, sd = 5e-5
        result = result * (dgamma(meanD,shape=0.02/0.018,scale=0.018) * normalization)
        }
    # penalty on meanO
    if(!missing(meanO) & !is.null(meanO))
	    {
	    # mean = 0.01, sd = 5e-5
	    result = result * (dgamma(meanO,shape=0.02/0.018,scale=0.018) * normalization)
	    }
  return(result)
}










    
