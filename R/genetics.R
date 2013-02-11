
known.alleles <- function(allNames, refData) {
  # Filters reference profile data to those used in this run.
  #
  # Also transforms data in profile to a vector of characters, e.g. "14,16" to
  # c("14", "16")
  #
  # Parameters:
  #   allNames: Name(s) of individuals we want in profile. The names should
  #   correspond to rows in refData. Although not specifically required by this
  #   function, the queried individual is expected to come first.
  # Returns: 
  #   reformats each locus by putting all alleles into a single vector.
  reformat.locus <- function(n) {
    result <- lapply(n, function(v) strsplit(v, ','))
    return(unlist(result))
  }
  # Now apply reformat.locus to each locus
  result <- lapply(refData[allNames, ], reformat.locus) 
  return(data.frame(result, stringsAsFactors=FALSE))
}

has.dropouts <- function(name, refData, cprofs) {
  # True if some alleles in profile are not CSP.
  #
  # Parameters:
  #   name: Name of the profile to check.
  #   refData: Reference profiles.
  #   cprofs: Crime Scene Profile, internal representation.
  queried = known.alleles(name, refData)[names(cprofs)]
  # Check dropouts for specific locus and replicate
  has.dropouts.per.rep = function(cprofsRep, queriedLocus) {
    validRep = !any(sapply(cprofsRep$csp, is.na))
    if(validRep && any(!queriedLocus %in% cprofsRep$csp)) return(TRUE)
    return(FALSE)
  }
  # Checks if there are any dropout for specific locus
  has.dropouts.per.locus = function(queriedLocus, cprofsLocus) {
    sapply(cprofsLocus, has.dropouts.per.rep, queriedLocus=queriedLocus) 
  }
  return( any(mapply(has.dropouts.per.locus, queried, cprofs) ))
}

known.without.dropouts <- function(name, refData, cprofs) {
  # Filters refDAta for those alleles without dropouts. 
  #
  # Returns in format of known.alleles
  dropouts = sapply(name, has.dropouts, refData=refData, cprofs=cprofs) 
  return(known.alleles(name[!dropouts], refData))
}
known.with.dropouts <- function(name, refData, cprofs) {
  # Filters refDAta for those alleles with dropouts. 
  #
  # Returns in format of known.alleles
  dropouts = sapply(name, has.dropouts, refData=refData, cprofs=cprofs) 
  return(known.alleles(name[dropouts], refData))
}
nb.with.dropouts <- function(nameK, refData, cprofs) {
  # Number of profiles with dropouts.
  #
  # Computes number of profiles named in nameK with dropouts. 
  #
  # Parameters:
  #   nameK: Names of the profile to check.
  #   refData: Reference profiles.
  #   cprofs: Crime Scene Profile, internal representation.
  return( sum(sapply( nameK, has.dropouts, refData=refData, cprofs=cprofs)) )
} 
nb.without.dropouts <- function(nameK, refData, cprofs) {
  # Number of profiles without dropouts.
  #
  # Computes number of profiles named in nameK with dropouts. 
  #
  # Parameters:
  #   nameK: Names of the profile to check.
  #   refData: Reference profiles.
  #   cprofs: Crime Scene Profile, internal representation.
  return(length(nameK) - nb.with.dropouts(nameK, refData, cprofs))
} 

ethnic.frequencies <- function(ethnic, loci=NULL, afreq=NULL) {
  # Reformats allele frequencies to include only given ethnic group and loci.
  #
  # Parameters:
  #   ethnic: Name of ethnic group. Should correspond to what's in afreq.
  #   loci: Names of the loci to use. Or use all.
  #   afreq: Frequency table. If NULL, loads in frequency table provided with
  #          likeLTD package.
 
  # Load frequency database if needed. 
  if(is.null(afreq)) afreq <- load.frequencies()
  # figueres out loci
  if(is.null(loci)) loci <- t(unique(afreq['Marker']))

  # Function which recreates the input for a single locus.
  # Input of a locus consists of one row per allele type.
  # Also filters out alleles with 0 frequencies 
  filter.locus <- function(n) {
    locus <- subset(afreq, Marker==n)
    result <- matrix(c(locus[[ethnic]], locus[["BP"]]), , 2)
    result[is.na(result[, 2]), 2] <- 0
    rownames(result) <- locus$Allele
    return(result[result[, 1] > 0, ])
  }
  # now apply function over all loci.
  result <- sapply(loci, filter.locus)
  # Then center around mean fragment length
  s1 = sum( sapply(result, function(n) sum(n[, 1] * n[, 2], na.rm=TRUE)) )
  s2 = sum( sapply(result, function(n) sum(n[, 1], na.rm=TRUE)) )
  for(j in 1:length(result)) result[[j]][,2] <- result[[j]][, 2] - s1/s2
  return(result)
}

add.missing.alleles <- function(frequencies, cprofs, queried, profiled=NULL) {
  # Add missing alleles to frequencies.
  #
  # The crime scene profile may present alleles which are missing from the
  # frequency database. We only add those alleles which present in queried
  # profile, and those in the known profiles which are subject to dropout.
  #
  # Parameters:
  #   frequencies: Arranged frequency table as returned by ethnic.frequencies.
  #   cprofs: Internal representation of the CSP.
  #   queried: queried profile as returned by known.alleles
  #   profiled: A list of loci containing further alleles. The CSP alleles
  #             which are not in these profiles and not in the database will be
  #             added to the database.
  # Returns: Frequency database filled with missing alleles. Also, it is
  #          reordered so htat loci follow the same order as cprofs.[
  for(locus in names(frequencies)) {
     freq.alleles <- row.names(frequencies[[locus]])

     # Alleles in queried profile not in database
     missing.alleles <- setdiff(queried[[locus]], freq.alleles)
     
     # Now check CSP alleles which are not in profiles and not in frequency
     # table.
     # Alleles which are in CSP. 
     csp.alleles <- unlist(sapply(cprofs[[locus]], function(n) n$csp))
     # Check that this is not an NA, 
     # And remove empty string.
     invalidCSP = unlist(sapply(cprofs[[locus]], function(n) is.na(n$csp)))
     if(any(invalidCSP)) csp.alleles <- c()
     else csp.alleles <- setdiff(csp.alleles, c(""))
     # skip if no alleles in CSP. 
     if(!setequal(csp.alleles, c())) {
       # Alleles which are in CSP but not in profiled
       csp.profiled <- setdiff(csp.alleles, profiled[[locus]])
       
       # Add to unknown alleles those in csp.profiled which are not in database
       missing.alleles <- union( missing.alleles, 
                                 setdiff(csp.profiled, freq.alleles) )
     }

     # Now add to database all missing alleles, with rownames
     if(!setequal(missing.alleles, c())) {
       names <- c(row.names(frequencies[[locus]]), missing.alleles)
       new.rows <- t(array(c(1, 0), c(2, length(missing.alleles))))
       frequencies[[locus]] <- rbind(frequencies[[locus]], new.rows)
       row.names(frequencies[[locus]]) <- names
     }
  } # loop over loci
  return(frequencies[names(cprofs)])
}

presence.matrices <- function(frequencies, cprofs, ..., type="csp") {
  # Matrices signifiying presence or absence of alleles in CSP.
  # 
  # A presence matrix for a given locus is an n by m matrix, where n is the
  # number of replicates and m is the number of alleles in the frequency
  # database. Each element (n,m) of the matrix is either TRUE if the allele is
  # present in that replicate of the CSP, or FALSE otherwise. 
  #
  # Note:
  #   Alleles which are in the CSP but not in the frequency database are
  #   *ignored*. It might be a good idea to call add.missing.alleles first.
  #
  # Note: 
  #  The following is expected to be TRUE
  #
  #  > names(frequencies) == names(cprofs) 
  #
  #  It has implications on the loci which are in the database and the order in
  #  which they are presented. 
  #
  # Parameters:
  #   frequencies: Allele frequency database as returned by ethnic.frequencies. 
  #   cprofs: CSP, internal representation
  #   ...:  Other lists over locus where each element is a list of
  #         alleles for which presence should be counted, much as is in cprofs.
  #         The outer (column) dimension of the list should be over loci, in
  #         the same order as in cprofs. The rows can be over replicates or
  #         something. In practice, these arguments are for non-dropout
  #         reference profiles in the presence matrix for "uncertain" alleles.
  #   type: Should be "csp" or "unc", depending on whether the presence matrix
  #         is computed for detected or uncertain alleles in CSP.
  # Returns: A list of presence matrices, one for each locus in CSP.
  # 
  # Note: The loci must occur in the same order in *all* input lists.  
  presence.per.locus <- function(frqLocus, cprofsLocus, ..., type="csp") {
    # Same as presence.matrices but for a single locus.
    
    # Check if valid locus first.
    invalidCSP = unlist(sapply(cprofsLocus, function(n) is.na(n$csp)))
    if( any(invalidCSP) ) return(NULL)

    # count: creates a row of the matrix.
    #        A row consists of 0 and 1, where 0 indicates that the allele is
    #        not present in the CSP. It is also possible to add other list of
    #        alleles via the ellipsis. In that case, a row consists of integers
    #        indicating the number of lists in which the allele is found. 
    rep = row.names(frqLocus)
    if(length(list(...)) == 0) {
      count <- function(n) as.integer(rep %in% n[[type]])
    }
    else {
      count <- function(n) {
        result <- c( as.integer(rep %in% n[[type]]),
                     sapply( list(...), function(n) as.integer(rep %in% n)) ) 
        result <- t(matrix(result, nrow=length(rep)))
        return(colSums(result))
      }
    }
    # create result matrix“
    result <- sapply(cprofsLocus, count) 
    result <- t(matrix(result, ncol=length(cprofsLocus)))
    return(result)
  }
  # Sanity check. Don't know how to do this for ellipsis though...
  if( any(names(frequencies) != names(cprofs)) ) {
    stop("CSP and frequencies do not have same names.")
  }

  return(mapply(presence.per.locus, frequencies, cprofs, ..., type=type))
}


adjust.frequencies <- function(frequencies, queriedAlleles, adj=1, fst=0.02) {
  # Adjust frequencies for current case.
  #
  # Parameters:
  #   frequencies: List of allele frequency per locus, for a given ethnic group
  #                (see ethnic.frequencies), completed for rare alleles in
  #                queried profile (see add.missing.alleles). The loci in both
  #                frequencies and queriedAlleles should occur in the same
  #                order.
  #   queriedAlleles: Profile of the queried individual, with colums as locus.
  #                   This will generally be something returned by
  #                   known.alleles.
  #   adj: sampling adjustment applied to the alleles of Q to avoid very low
  #        counts for rare alleles, and to allow for the allele counts to take
  #        into account the genotype of Q.  The default value is adj = 1; there
  #        are reasons to prefer adj = 2 but the value of adj is much less
  #        important than Fst (see below) unless the database size is very
  #        small.  
  #   fst: allows for shared ancestry of Q with X. Recommended that Fst should
  #        be at least 0.02, and may need to be as high as 0.05 in some
  #        populations (e.g.  small, isolated subpopulations of the population
  #        from which the reference database has been drawn).

  adjust.per.locus <- function(frqLocus, q.loc, adj=1, fst=0.02) {
    # Applies operations to a single locus.
    homozygote <- as.integer(q.loc[1] == q.loc[2])
    frqLocus[q.loc, 1] <- frqLocus[q.loc, 1] + adj * (1 + homozygote)
    #  shared ancestry adjustements.
    frqLocus[, 1] <- frqLocus[, 1] / sum(frqLocus[, 1]) * (1 - fst) / (1 + fst) 
    frqLocus[q.loc, 1] <- frqLocus[q.loc, 1] + fst / (1 + fst) * (1 + homozygote)
    return(frqLocus)
  }
  return(mapply(adjust.per.locus, frequencies, queriedAlleles, adj=adj,
                fst=fst))
}

all.profiles.per.locus = function(nAlleles, nContrib=1) {
  # All possible agreggate profiles for a given locus.
  #
  # Computes all possible combined profiles from n contributors.
  # This is a permutation on combinations: it grows way too fast as
  # :math:`N=[0.5*nAlleles*(nAlleles+1)]^{nContrib}`.
  # 
  # Parameters:
  #   nAlleles: number of alleles. 
  #   nContrib: number of unknowns contributors.
  
  # All possible genotypes for a single contributor.
  singleContributor = combinations(nAlleles, 2, rep=T)
  if(nContrib == 1) return(singleContributor)
  
  # All possible permutations of "nContrib" contributors.
  nContribPerms = permutations(nrow(singleContributor), nContrib, rep=T)
  # Function to apply single column from contributor permutations to the allele
  # combinations of a single contributor.
  apply.perms <- function(n, alleles, perms) alleles[perms[, n], ]
  # All possible allele combintations of "nContrib" contributors.
  unknownContributors <- do.call( cbind, 
                                  lapply( 1:nContrib, apply.perms, 
                                          alleles=singleContributor, 
                                          perms=nContribPerms ) )
  return(unknownContributors)
}

genoperms.per.locus = function(cspPresence, profPresence,
                               alleleNames, nUnknowns, dropin) {
  # Profiles of unknown contributors for given locus.
  #
  # Figures out all possible profiles of n contributors, including dependance
  # on doprin.
  #
  # Parameters:
  #   cspPresence: Crime Scene Profile for given locus. This should be a matrix
  #             where rows are replicas and contributors, columns are alleles
  #             for the given locus, and the input is TRUE or FALSE and
  #             indicates the presence of that particular allele.
  #   profPresence: Known profiles for given locus. Same type of format as
  #              cspLocus.
  #   alleleNames: Name of the alleles, e.g. columns of the previous two.
  #   nUnknowns: number of unknown contributors. 
  #   dropin: TRUE if modelling drop-ins.
  #
  # Return: A m by (2n) matrix where the colums (grouped by twos) correspond to
  #         contributors, and each row is their potential contribution to the
  #         CSP.
  #
  # Note: If nUnknowns == 0 and dropin == FALSE, then returns a matrix with the
  #       alleles in the queried profile. This is a logic different from the
  #       main case. 
  #
  # Note: If nUnknowns == 0 and dropin == TRUE is equivalent to nUnknowns == 1
  #       and dropin == TRUE. 

  # Case where there are no unknown contributors and no dropin.
  # Then the result will be the genotype of queried only.
  if(nUnknowns == 0 && dropin == FALSE) {
    cspPresence = cspPresence[nrow(cspPresence)-1: nrow(cspPresence), ]
    if(is.matrix(cspPresence)) cspPresence <- colSums(cspPresence) > 0
    else                       cspPresence <- cspPresence > 0

	  return(t(matrix(which(cspPresence>0))))
  }
  
  # Case where there are no unknown contributors but dropin is modelled.
  # Seems it's formally equivalent to one unknown and dropin.
  if(nUnknowns == 0 && dropin == TRUE) nUnknowns = 1

  # Case without dropin: check that there are enough unknown contributors to
  # account for the alleles in the CSP that are not in the known profile.
  # Might be able to return early here also.
  if(dropin == FALSE) {
    # cspPresence: logical vector indicating alleles occurence in CSP.
    if(is.matrix(cspPresence)) cspPresence <- colSums(cspPresence) > 0
    else                       cspPresence <- cspPresence > 0
    
    # profPresence: logical vector indicating alleles occurence in aggregate
    # of all input profile
    if(is.matrix(profPresence)) profPresence <- colSums(profPresence) > 0
    else                        profPresence <- profPresence > 0
    
    # knownUnknowns: indices of the alleles in the crime scene which are not
    #                in the known profiles. Indices are for freqLocus rows.
    #                in the case of dropins, then there are no known unknowns.
    knownUnknowns = which(cspPresence & !profPresence) 

    # Not enough contributors, return empty matrix.
    if(length(knownUnknowns) > 2*nUnknowns)  
      stop("Not enough unknown contributors.")
    if(length(knownUnknowns) == 2*nUnknowns) return(matrix(0, 1, 0))
  }

  # Compute all genotype permutations. 
  unknownContributors = all.profiles.per.locus(length(alleleNames), nUnknowns)

  # If dropins are not modelled, then make sure that alleles in CSP but not in
  # profiled are in (sum of) unknown contributors. 
  if (dropin == FALSE) {
    hasKnownUnknowns <- apply( unknownContributors, 1, 
                               function(n) all(knownUnknowns %in% n) )
    unknownContributors <- unknownContributors[hasKnownUnknowns, , drop=FALSE]
  }

  if(is.matrix(unknownContributors) == FALSE)
    unknownContributors <- t(as.matrix(unknownContributors))

  return(unknownContributors) 
}

#-----------------------------------------
# calc.fixed()
#-----------------------------------------
# Calculates all objects that vary by locus, but fixed by iteration. Differs for pros and def.

# Arguments required are:
#	Nkdo: Number of knowns subject to dropout. Set in global.objects()
#	Nunp: number of unprofiled contributors, either NU (pros) or NU+1 (def). Set in GUI
#	afbp: Adjusted allele fractions and allele lengths. Use af from prepro()
#	known.alleles: Use known[[j]][mfunpr]. known and mfunpr from global.objects()
#	CSP: CSP for all replicates. Use csp from prepro(). Different for pros and def if Q is not subject to dropout.
#	DI: DI probability (used to determine whether to use Ureq or not), from Drin. Set in GUI.
#	Qcont:  Use 0 for def, or Qdrop or 1+Qdrop for pros. Set in global.objects(). 0 if Q is not assumed to be a contributor (i.e. under Hd), =1 if Q is a contributor not subject to dropout and =2 if Q is a contributor subject to dropout.

# Returns:
#	kpdo: Alleles form known contributors. Differ for pros and def ans def cannot include Q's alleles. Also specific order set.
#	pUall: Generates all the permutations of genotypes for the unprofiled contributors.
#	index: internal abstraction
#	fragments: matrix of fragment lengths across all permutations
#	v.index: Internal abstraction required to index specific elements in tprof matrix.

calc.fixed = function(Nkdo,Nunp,afbp,known.alleles,CSP,DI,Qcont){

	# generates kpdo
	Qgen=known.alleles[1:2]
	if(Qcont==0)kpdo=known.alleles[-(1:2)] # defence case. Excludes Qs alleles.
	if(Qcont==1)kpdo=known.alleles # prosecution case. Includes Qs alleles.
	if(Qcont==2)kpdo=c(known.alleles[-(1:2)],Qgen) # if Q is contributor subject to dropout, put those alleles last in kpdo 

	# generates presence/absence of known alleles
	presence = rep(0,nrow(afbp)) 
	if(length(kpdo)) for(u in 1:length(kpdo))presence = presence + (rownames(afbp)==kpdo[u])

if(Nunp>0) {
	# generates pUall
	if(nrep>1) CSPset = colSums(CSP[CSP[,1]<999,])>0  else CSPset = CSP # T for alleles that occur in CSP at least once

	if (DI==0) {
		Ureq = which((CSPset * (presence == 0)) > 0)  # identify alleles in CSP but not in kpdo, so must come from the U or X
		nUx =  2*Nunp - length(Ureq)
  } else {nUx <- 2*Nunp}
  if(nUx < 0) {return(0)} else if(nUx > 0) {Ux = combinations(nrow(afbp),nUx,rep=T)} else {Ux = matrix(0,1,0)} 
		pUall.end <- matrix(data=c(rep(0, times=2*Nunp)),ncol=2*Nunp)
		for(j in 1:nrow(Ux)){  # all possible allocations of required + random alleles to the unprofileds
			if (DI==0){
        pUall =
				if(length(unique(c(Ureq,Ux[j,])))==(2*Nunp)) {permutations(2*Nunp,2*Nunp,c(Ureq,Ux[j,]),set=F)} else {unique(permutations(2*Nunp,2*Nunp,c(Ureq,Ux[j,]),set=F))}
				}
			else {pUall = 
				if(length(unique(Ux[j,]))==(2*Nunp)) {permutations(2*Nunp,2*Nunp,Ux[j,],set=F)}else {unique(permutations(2*Nunp,2*Nunp,Ux[j,],set=F))} 
				}
			
			pUall = pUall[which(apply(pUall[,2*(1:Nunp)-1,drop=F]<=pUall[,2*(1:Nunp),drop=F],1,prod)>0),,drop=F]
			pUall.end <- rbind(pUall.end,pUall)
			}
	pUall <- pUall.end[-1,]
} 


if(Nunp==0) {
	if(DI==0) {
		ind=c(1,1)
		if(nrep>1) CSPset = colSums(CSP[CSP[,1]<999,])>0  else CSPset = CSP
		presenceQ = rep(0,nrow(afbp)) 
		for(u in 0:1) {presenceQ = presenceQ + (rownames(afbp)==kpdo[length(kpdo)-u])}
		pUall <- which(presenceQ>0)   # If no unknowns set pUall to genotype of Q
		} 

	if(DI!=0) {
		ind=c(1,1)
		nUx <- 2
		if(nUx < 0) {return(0)} else if(nUx > 0) {Ux = combinations(nrow(afbp),nUx,rep=T)} else {Ux = matrix(0,1,0)} 
		pUall.end <- matrix(data=c(rep(0, times=length(ind))),ncol=length(ind))
		for(j in 1:nrow(Ux)){ 
			pUall = if(length(unique(Ux[j,]))==(length(ind))) {permutations(length(ind),length(ind),Ux[j,],set=F)}else {unique(permutations(length(ind),length(ind),Ux[j,],set=F))} 
			pUall = pUall[which(apply(pUall[,1,drop=F]<=pUall[,2,drop=F],1,prod)>0),,drop=F]
			pUall.end <- rbind(pUall.end,pUall)
			}
		pUall <- pUall.end[-1,]
		}
	}


# CHANGE: edit pUall
if(is.matrix(pUall)==FALSE) pUall <- t(as.matrix(pUall))

	# generates other fixed values
	N = dim(pUall)[1] # number of genotypes 
	fragments = t(matrix(afbp[pUall,2],nrow=N))
	v.index = (0:(N-1))*nrow(afbp)
	if(Nunp>0) index = (length(kpdo)/2) + rep((1:Nunp),rep(2,Nunp))  else index = (length(kpdo)/2) + ind


return(list(kpdo=kpdo,pUall=pUall,index=index,fragments=fragments,v.index=v.index))}
#-----------------------------------------
