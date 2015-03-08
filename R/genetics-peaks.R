allExplained = function(genotype,cspAlleles,knownWithStutter,alleleNames,doDoubleStutter=FALSE,doOverStutter=FALSE)
	{
	genotypeAlleles = as.numeric(alleleNames[genotype])
	genWithStutter = c(genotypeAlleles,genotypeAlleles-1)
        if(doDoubleStutter) genWithStutter = c(genWithStutter,genotypeAlleles-2)
	if(doOverStutter) genWithStutter = c(genWithStutter,genotypeAlleles+1)
	all(cspAlleles%in%c(genWithStutter,knownWithStutter))
	}

explain.all.peaks = function(cspPresence,profPresence,knownProfs,alleleNames,nUnknowns,cspAlleles,cspHeights, doDoubleStutter=FALSE,doOverStutter=FALSE)
	{
  	# Catch early. Shouldn't be a problem here, but will be at some point. 
	if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
	if(!is.matrix(profPresence)) stop("Expected a matrix as input.")
	if(ncol(profPresence)!=0)
	    {
	    # get known alleles
    	knownIndex = sapply(unlist(knownProfs),FUN=function(x) which(alleleNames==x))
    	knownAlleles = as.numeric(alleleNames[knownIndex])
	# include stuttered known alleles
	knownWithStutter = c(knownAlleles,knownAlleles-1)
	if(doDoubleStutter) knownWithStutter = c(knownWithStutter, knownAlleles-2)
	if(doOverStutter) knownWithStutter = c(knownWithStutter, knownAlleles+1)
	    } else {
        knownIndex = c()
        knownWithStutter = c()
	    }
	# get csp Alleles
	cspAlleles = alleleNames[row(cspPresence)[which(cspPresence)]]
	# if no unknowns return knowns
	if(nUnknowns==0) 
		{
		check = !all(round(as.numeric(cspAlleles),1)%in%round(as.numeric(knownWithStutter),1))
		if(check) 
		    {
		    stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
		    } else {
		    return(matrix(knownIndex,ncol=1))
		    }
		}
	# get all genotype combinations for unknowns
	genCombs = all.genotypes.per.locus(length(alleleNames),nUnknowns)
	# find which combinations explain all peaks
	index = apply(genCombs,MARGIN=2,FUN=function(x) allExplained(x,cspAlleles,knownWithStutter,alleleNames,doDoubleStutter,doOverStutter))
	if(length(which(index))==0) stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
	genCombs = genCombs[,index,drop=FALSE]

	# add known profiles as last contributors
	if(length(knownIndex)!=0) 
		{
		genCombs = rbind(genCombs,matrix(rep(knownIndex,times=ncol(genCombs)),ncol=ncol(genCombs)))
		}
	# remove very unlikely genotype combinations
	#index = apply(matrix(as.numeric(alleleNames[genCombs]),ncol=ncol(genCombs)),MARGIN=2,FUN=function(x) any(likeLTD:::unlikelyGenotypes(x,cspAlleles,cspHeights)))
	#print(paste0(length(which(index))," removed out of ", ncol(genCombs)))
	#if(length(which(index))!=0&length(which(index))!=ncol(genCombs)) genCombs = genCombs[,index,drop=FALSE]
	return(genCombs)
	}

unlikelyGenotypes = function(genotype,cspAlleles,cspHeights)
	{
	# gen = single allele in genotype, genotype = total genotype
	checkRemove = function(gen,genotype,cspAlleles,cspHeights)
		{
		if(gen==-1) return(FALSE)
		if(!gen%in%cspAlleles&(gen-1)%in%cspAlleles&!(gen-1)%in%genotype) return(TRUE) 
		# taking into account double stutter
		if(!gen%in%cspAlleles&(gen-2)%in%cspAlleles&!(gen-1)%in%cspAlleles&!(gen-2)%in%genotype) return(TRUE) else return(FALSE)
		}
	out = sapply(genotype,FUN=function(x) checkRemove(x,genotype,cspAlleles,cspHeights))
	return(out)
	}

# Profiles of unknown contributors for given locus.
compatible.genotypes.peaks = function(cspPresence, profPresence, knownProfs,alleleNames,
                                nUnknowns, dropin=FALSE, missingReps=NULL) {
  
  # Catch early. Shouldn't be a problem here, but will be at some point. 
  if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
  if(!is.matrix(profPresence)) stop("Expected a matrix as input.")

  # Case where there are no unknown contributors but dropin is modelled.
  # Seems it's formally equivalent to one unknown and dropin.
  if(nUnknowns == 0 && dropin == TRUE) nUnknowns = 1

  # Compute all genotype permutations. 
  genotypes = all.genotypes.per.locus(length(alleleNames), nUnknowns)
  if(dropin) return(genotypes) # Dropin case: can return immediately.

  # Case without dropin: check that there are enough unknown contributors to
  # account for the alleles in the CSP that are not in the known profile.
  # Might be able to return early here also.
  # Reduce matrices to a single logical vector
  if(is.null(missingReps)) missingReps = rep(TRUE, nrow(cspPresence))
  cspPresence  = rowSums(cspPresence[, c(!missingReps), drop=FALSE]) > 0
  profPresence = rowSums(profPresence) > 0

  # required: indices of the alleles in the crime scene which are not in the
  #           known profiles. Indices are for freqLocus rows. In the case of
  #           dropins, then there are no required alleles.
  required = which(cspPresence & !profPresence) 

  # Not enough contributors, return empty matrix.
  if(length(required) > 2*nUnknowns)  {
    stop(sprintf("Not enough unknown contributors:
       %d contributors, 
       %d required allele (%s).", 
       nUnknowns, length(required),
       paste(alleleNames[required], collapse=", ")))
  }
  #if(nUnknowns == 0) return(matrix(0, nrow=1, ncol=1))


  hasRequired <- apply(genotypes, 2, function(n) all(required %in% n))
  genotypes = genotypes[, hasRequired, drop=FALSE]

# if no unknown contributors, sometimes have an empty genotypes object
if(ncol(genotypes)==0) genotypes = matrix(ncol=1,nrow=0)

  # add known profiles as last contributors
  knownAlleles = sapply(unlist(knownProfs),FUN=function(x) which(alleleNames==x))
  if(length(knownAlleles)!=0) genotypes = rbind(genotypes,matrix(rep(knownAlleles,times=ncol(genotypes)),ncol=ncol(genotypes)))

return(genotypes) 
}

# Defines prod.matrix from R.
prod.matrix.row <- function(x) {
  # Fast column-wise product.
  #
  # Sadly, this is faster than apply(x, 1, prod)
  y=x[,1]
  for(i in 2:ncol(x))
  y=y*x[, i]
  return(y)
}

ethnic.database.lus <- function(ethnic, loci=NULL, afreq=NULL) {
  # Reformats allele database to include only given ethnic group and loci.
  #
  # Filters database down to a single ethnic group and the loci of interests.
  # Removes alleles which are not present in given ethnic group. Centers
  # average length across filtered database.
  # 
  # Parameters:
  #   ethnic: Name of ethnic group. Should correspond to what's in afreq.
  #   loci: Names of the loci to use. Or use all.
  #   afreq: Frequency table. If NULL, loads in frequency table provided with
  #          likeLTD package.
 
  # Load frequency database if needed. 
  if(is.null(afreq)) afreq <- load.allele.database()
  # figueres out loci
  if(is.null(loci)) loci <- t(unique(afreq['Marker']))

  # Function which recreates the input for a single locus.
  # Input of a locus consists of one row per allele type.
  # Also filters out alleles with 0 frequencies.
  filter.locus <- function(n) {
    locus <- afreq[afreq$Marker == n, ]
    result <- matrix(c(locus[[ethnic]], locus[["BP"]], locus[["LUS"]]), , 3)
    result[is.na(result[, 2]), 2] <- 0
    rownames(result) <- locus$Allele
    result = fill.unknown.LUS(result)
    return(result[result[, 1] > 0, ])
    return(result)
  }
  # now apply function over all loci.
  result <- sapply(loci, filter.locus)
  # Then center around mean fragment length
  s1 = sum( sapply(result, function(n) sum(n[, 1] * n[, 2], na.rm=TRUE)) )
  s2 = sum( sapply(result, function(n) sum(n[, 1], na.rm=TRUE)) )
  for(j in 1:length(result)) result[[j]][, 2] <- result[[j]][, 2] - s1/s2
  return(result)
}
