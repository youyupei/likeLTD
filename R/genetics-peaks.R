allExplained = function(genotype,cspAlleles,knownWithStutter,alleleNames)
	{
	genotypeAlleles = as.numeric(alleleNames[genotype])
genWithStutter = unique(c(genotypeAlleles,genotypeAlleles-1))
	all(cspAlleles%in%c(genWithStutter,knownWithStutter))
	}

explain.all.peaks = function(cspPresence,profPresence,knownProfs,alleleNames,nUnknowns,cspAlleles,cspHeights)
	{
	CSPPRESENCE <<- cspPresence
PROFPRESENCE <<- profPresence
KNOWNPROFS <<- knownProfs
ALLELENAMES <<- alleleNames
NUNKNOWNS <<- nUnknowns
CSPALLELES <<- cspAlleles
CSPHEIGHTS <<- cspHeights
  	# Catch early. Shouldn't be a problem here, but will be at some point. 
	if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
	if(!is.matrix(profPresence)) stop("Expected a matrix as input.")
	if(ncol(profPresence)!=0)
	    {
	    # get known alleles
    	knownIndex = sapply(unlist(knownProfs),FUN=function(x) which(alleleNames==x))
    	knownAlleles = as.numeric(alleleNames[knownIndex])
    	# include stuttered known alleles
    	knownWithStutter = unique(c(knownAlleles,knownAlleles-1))
	    } else {
        knownIndex = c()
        knownWithStutter = c()
	    }
	# get csp Alleles
	cspAlleles = alleleNames[row(cspPresence)[which(cspPresence)]]
	# if no unknowns return knowns
	if(nUnknowns==0) 
		{
		if(!all(cspAlleles%in%knownWithStutter)) stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
		return(matrix(knownIndex,ncol=1))
		}
	# get all genotype combinations for unknowns
	genCombs = likeLTD:::all.genotypes.per.locus(length(alleleNames),nUnknowns)
	# find which combinations explain all peaks
	index = apply(genCombs,MARGIN=2,FUN=function(x) likeLTD:::allExplained(x,cspAlleles,knownWithStutter,alleleNames))
	if(length(which(index))==0) stop(paste0("Not enough contributors to explain CSP at locus ", colnames(knownProfs)))
	genCombs = genCombs[,index]

	# add known profiles as last contributors
	if(length(knownIndex)!=0) 
		{
		genCombs = rbind(genCombs,matrix(rep(knownIndex,times=ncol(genCombs)),ncol=ncol(genCombs)))
		}
	# remove very unlikely genotype combinations
	index = apply(matrix(as.numeric(alleleNames[genCombs]),ncol=ncol(genCombs)),MARGIN=2,FUN=function(x) any(likeLTD:::unlikelyGenotypes(x,cspAlleles,cspHeights)))
	#print(paste0(length(which(index))," removed out of ", ncol(genCombs)))
	if(length(which(index))!=0&length(which(index))!=ncol(genCombs)) genCombs = genCombs[,index]
	return(genCombs)
	}

unlikelyGenotypes = function(genotype,cspAlleles,cspHeights)
	{
	checkRemove = function(gen,genotype,cspAlleles,cspHeights)
		{
		if(gen==-1) return(FALSE)
		if(!gen%in%cspAlleles&(gen-1)%in%cspAlleles&!(gen-1)%in%genotype) return(TRUE) else return(FALSE)
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
  genotypes = likeLTD:::all.genotypes.per.locus(length(alleleNames), nUnknowns)
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
