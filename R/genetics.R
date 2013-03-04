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

ethnic.database <- function(ethnic, loci=NULL, afreq=NULL) {
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
  # Also filters out alleles with 0 frequencies 
  filter.locus <- function(n) {
    locus <- afreq[afreq$Marker == n, ]
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
  for(j in 1:length(result)) result[[j]][, 2] <- result[[j]][, 2] - s1/s2
  return(result)
}

add.missing.alleles <- function(alleleDb, cprofs, queried, profiled=NULL) {
  # Add missing alleles to alleleDb
  #
  # The crime scene profile may present alleles which are missing from the
  # frequency database. We only add those alleles which present in queried
  # profile, and those in the known profiles which are subject to dropout.
  #
  # Parameters:
  #   alleleDb: Arranged frequency table as returned by ethnic.database
  #   cprofs: Internal representation of the CSP.
  #   queried: queried profile as returned by known.alleles
  #   profiled: A list of loci containing further alleles. The CSP alleles
  #             which are not in these profiles and not in the database will be
  #             added to the database.
  # Returns: Frequency database filled with missing alleles. Also, it is
  #          reordered so htat loci follow the same order as cprofs.
  for(locus in names(alleleDb)) {
     freq.alleles <- row.names(alleleDb[[locus]])

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
       # Alleles which are in CSP but not in profiled
       csp.profiled <- setdiff(csp.alleles, profiled[[locus]])
       
       # Add to unprofiled alleles those in csp.profiled which are not in database
       missing.alleles <- union( missing.alleles, 
                                 setdiff(csp.profiled, freq.alleles) )
     }

     # Now add to database all missing alleles, with rownames
     if(!setequal(missing.alleles, c())) {
       names <- c(row.names(alleleDb[[locus]]), missing.alleles)
       new.rows <- t(array(c(1, 0), c(2, length(missing.alleles))))
       alleleDb[[locus]] <- rbind(alleleDb[[locus]], new.rows)
       row.names(alleleDb[[locus]]) <- names
     }
  } # loop over loci
  return(alleleDb[names(cprofs)])
}

missing.csp.per.locus <- function(cprofsLocus) {
  # Logical array indicating missing replicates with TRUE.
  #
  # Parameters:
  #   cprofsLocus: cprofs for a given locus.
  return(sapply(cprofsLocus, function(n) length(n) <= 1)) 
}

  
presence.matrices <- function(alleleDb, cprofs, ..., type="csp") {
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
  #  > names(alleleDb) == names(cprofs) 
  #
  #  It has implications on the loci which are in the database and the order in
  #  which they are presented. 
  #
  # Parameters:
  #   alleleDb: Allele frequency database as returned by ethnic.database. 
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
  presence.per.locus <- function(alleleDbLocus, cprofsLocus, ..., type="csp") {
    # Same as presence.matrices but for a single locus.
    
    # Check if valid locus first.
    if( any(missing.csp.per.locus(cprofsLocus)) ) return(NULL)

    # count: creates a row of the matrix.
    #        A row consists of 0 and 1, where 0 indicates that the allele is
    #        not present in the CSP. It is also possible to add other list of
    #        alleles via the ellipsis. In that case, a row consists of integers
    #        indicating the number of lists in which the allele is found. 
    rep = row.names(alleleDbLocus)
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
  if( any(names(alleleDb) != names(cprofs)) ) {
    stop("CSP and alleleDb do not have same names.")
  }

  return(mapply(presence.per.locus, alleleDb, cprofs, ..., type=type))
}


adjust.frequencies <- function(alleleDb, queriedAlleles, adj=1, fst=0.02) {
  # Adjust frequencies for current case.
  #
  # Parameters:
  #   alleleDb: Allele database for given locus, for a given ethnic group (see
  #             ethnic.database), completed for rare alleles in queried
  #             profile (see add.missing.alleles). The loci in both alleleDb
  #             and queriedAlleles should occur in the same order.
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
  adjust.per.locus <- function(alleleDbLocus, queriedLocus, adj=1, fst=0.02) {

    # Applies operations to a single locus.
    homozygote <- as.integer(queriedLocus[1] == queriedLocus[2])
    alleleDbLocus[queriedLocus, 1] <- alleleDbLocus[queriedLocus, 1] + adj * (1 + homozygote)
    #  shared ancestry adjustements.
    alleleDbLocus[, 1] <-
      alleleDbLocus[, 1] / sum(alleleDbLocus[, 1]) * (1 - fst) / (1 + fst) 
    alleleDbLocus[queriedLocus, 1] <-
      alleleDbLocus[queriedLocus, 1] + fst / (1 + fst) * (1 + homozygote)
    return(alleleDbLocus)
  }
  return(mapply(adjust.per.locus, alleleDb, queriedAlleles, adj=adj,
                fst=fst))
}

all.genotypes.per.locus = function(nAlleles, nContrib=1) {
  # All possible agreggate profiles for a given locus.
  #
  # Computes all possible combined profiles from n contributors.
  # This is a permutation on combinations: it grows way too fast as
  # :math:`N=[0.5*nAlleles*(nAlleles+1)]^{nContrib}`.
  # 
  # Parameters:
  #   nAlleles: number of alleles. 
  #   nContrib: number of unprofiled contributors.
  
  # All possible genotypes for a single contributor.
  singleContributor = t(combinations(nAlleles, 2, repeats.allowed=T))
  if(nContrib == 1) return(singleContributor)
  
  # All possible permutations of "nContrib" contributors.
  nContribPerms = permutations(ncol(singleContributor), nContrib,
                               repeats.allowed=T)
  # The next two lines create a matrix where the first two columns is
  # singleContributor in the order given by the first column of nContribPerms.
  # The next two columns is singleContributor in the order given by
  # nContribPerms's second column. And so on and so forth.
  apply.perms <- function(n, sing, perms) sing[, perms[, n]]
  rows <- lapply(1:nContrib, apply.perms, sing=singleContributor,
                 perms=nContribPerms) 
  do.call(rbind, rows)
}

possible.genotypes = function(cspPresence, profPresence, missingReps,
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
  # Note: If nUnknowns == 0 and dropin == TRUE is equivalent to nUnknowns == 1
  #       and dropin == TRUE. 

  
  # Catch early. Shouldn't be a problem here, but will be at some point. 
  if(!is.matrix(cspPresence)) stop("Expected a matrix as input.")
  if(!is.matrix(profPresence)) stop("Expected a matrix as input.")

  # Case where there are no unknown contributors but dropin is modelled.
  # Seems it's formally equivalent to one unknown and dropin.
  if(nUnknowns == 0 && dropin == TRUE) nUnknowns = 1

  # Compute all genotype permutations. 
  genotypes = all.genotypes.per.locus(length(alleleNames), nUnknowns)
  if(dropin) return(genotypes) # Dropin case: can return immediately.

  # Case without dropin: check that there are enough unknown contributors to
  # account for the alleles in the CSP that are not in the known profile.
  # Might be able to return early here also.
  # Reduce matrices to a single logical vector
  cspPresence  = rowSums(cspPresence[, c(!missingReps), drop=FALSE]) > 0
  profPresence = rowSums(profPresence) > 0

  # required: indices of the alleles in the crime scene which are not in the
  #           known profiles. Indices are for freqLocus rows. In the case of
  #           dropins, then there are no required alleles.
  required = which(cspPresence & !profPresence) 

  # Not enough contributors, return empty matrix.
  if(length(required) > 2*nUnknowns)  stop("Not enough unknown contributors.")
  if(nUnknowns == 0) return(matrix(0, nrow=1, ncol=1))
  if(length(required) == 2*nUnknowns) return(matrix(0, 1, 0))

  hasRequired <- apply(genotypes, 2, function(n) all(required %in% n))
  genotypes[, hasRequired, drop=FALSE]
}

relatedness <- function(nAlleles, nContribs, profiles, ibds) {
  # Compute relatedness profiles
  # 
  # There are three profiles: one where the first allele is the first IBD, one
  # wher the first alllele is always the second IBD, one where a whole
  # contributor must be the IBD.

  relprofile <- list()
  for(i in 1:2) {
    condition <- apply(profiles, 2, function(n) ibds[i] %in% n[1:2])
    relprofile[[i]] <- profiles[, condition]
    condition <- relprofile[[i]][2, ] == ibds[i]
    relprofile[[i]][2, condition] <- relprofile[[i]][1, condition]
    relprofile[[i]][1, condition] <- ibds[i]
  }
  
  relprofile[[3]] <- all.genotypes.per.locus(nAlleles, nContribs-1)
  relprofile[[3]] <- rbind(ibds[1], ibds[2], relprofile[[3]])
  relprofile
}

known.epg.per.locus <- function(relContrib, degradation, fragmentLengths,
                                profiles) {
  # Creates "electropherogram" vector for a given locus.
  #
  # In practice, the "electropherogram" is a vector giving the dose for each
  # allele present in the known profiles.
  #
  # Parameters:
  #   relContrib: relative contributions from each individual. It should be
  #               vector of the same length as the number of individuals in the
  #               profile.
  #   degradation: degradation of the DNA from each individual. 
  #   fragmentLengths: The str lengths for a given locus.
  #   profiles: A matrix of m by (2n), where is the number of individuals in
  #             the profile, m is the number of alleles in the frequency table,
  #             and each element indicate the presence or absence of that
  #             particular allele in the makeup of each individual.
  # Returns: A vector with an element per allele. Each element is zero if is
  #          not within the known profile, or it is its dose within the
  #          candidate CSP.

  indices = which(profiles > 0) - 1
  indivs = trunc(indices / nrow(profiles)) + 1
  alleles = indices %% nrow(profiles) + 1
  doses = mapply( function(d, r, L, het) het * r * (1.0 + d)^{-L}, 
                  degradation[indivs], 
                  relContrib[indivs], 
                  fragmentLengths[alleles],
                  rep(3 - colSums(profiles), colSums(profiles)) )
  result = array(0.0, nrow(profiles))
  for(i in 1:length(alleles)) 
    result[alleles[i]] = result[alleles[i]] + doses[i]
  return(result)
}

all.epg.per.locus <- function(relContrib, degradation, profPresence,
                              knownFragLengths, fragLengths, genotypes,
                              addProfiles) {
  # Creates "electropherogram" for each possible set of unknown contributors.
  #
  # Parameters:
  #   relContrib: relative contributions from each individual. It should be
  #               vector of the same length as the number of individuals in the
  #               profile.
  #   degradation: degradation of the DNA from each individual. 
  #   profPresence: Presence matrix for known profiles.
  #   knownFragLengths: matrix with fragment lengths for each allele in each
  #                     known profile.
  #   fragLengths: matrix with fragment lengths for each allele each computed
  #                profile.
  #   genotypes: matrix with all possible genotypes. 
  # Returns: A vector with an element per allele. Each element is zero if is
  #          not within the known profile, or it is its dose within the
  #          candidate CSP.

  # "electropherogram" from known profiles
  knownEPG <- known.epg.per.locus(relContrib, degradation, knownFragLengths,
                                  profPresence)

  # allEPG: for each set of unknown contributors, total EPG of known and
  #         unknowns.
  # At point of creation, contains only component from known profiles.
  allEPG = matrix(knownEPG, ncol=ncol(genotypes), nrow=length(knownEPG)) 
  # Add into allEPG the components from unknown contributors.
  if(addProfiles) {
    nUnknowns = nrow(genotypes)
    # v.index: index matrix to access alleles across possible agreggate
    # profiles. 
    v.index = 0:(ncol(genotypes)-1) * length(knownFragLengths)
    indices = ncol(profPresence) + rep(1:(nUnknowns/2), rep(2, nUnknowns/2))
    unknownDoses = relContrib[indices] *
                   (1.0 + degradation[indices])^-fragLengths
    # Perform sum over each DNA strand of each unknown contributor.
    for(u in 1:nUnknowns){
      indices = genotypes[u, ] + v.index
      allEPG[indices] = allEPG[indices]  + unknownDoses[u, ]
    }
  }
  return(allEPG)
}

# Defines prod.matrix from R.
prod.matrix.R <- function(x) {
  # Fast row-wise product
  #
  # Sadly, this is faster than apply(x, 1, prod)  
  y=x[,1]
  for(i in 2:dim(x)[2])
  y=y*x[,i]
  return(y)
}
# Tries and defines a fast prod.matrix if possible.
if(require("inline")) {
  code = "SEXP res;
          int const nrow = INTEGER(GET_DIM(x))[0];
          int const ncol = INTEGER(GET_DIM(x))[1];
          PROTECT(res = allocVector(REALSXP, nrow));
          double *inptr = REAL(x);
          double *const outptr_first = REAL(res);
          double *outptr = outptr_first;
          // Copy first column into output vector.
          for(int i=0; i < nrow; ++i, ++inptr, ++outptr)
            *outptr = (*inptr == NA_INTEGER) ? 1: *inptr;
          // Perform multiplication with other columns
          for(int j=1; j < ncol; ++j) {
            outptr = outptr_first;
            for(int i=0; i < nrow; ++i, ++inptr, ++outptr)
              if(*inptr != NA_REAL) *outptr *= *inptr;
          }
          UNPROTECT(1);
          return res;"
  prod.matrix.C = cfunction(signature(x="array"), code)
  prod.matrix = prod.matrix.C
} else { prod.matrix = prod.matrix.R }
prod.matrix = prod.matrix.R

selective.row.prod <- function(condition, input) {
  # Row-wise product over selected columns or elements.
  #
  # If all(condition == FALSE), returns 1. Otherwise does one of the following.
  #
  # - If condition and input matrices: 
  #   Performs column-wise multiplication of input elements, ignoring those
  #   which are *not* selected by the condition matrix.
  # - If condition is a matrix and input a vector: 
  #   Equivalent to creating a matrix whith the same dimensions as condition
  #   from the input, where row ``i`` is ``input[condition[i, ]]``, then apply
  #   product.
  # - If condition is a vector and input is a matrix: 
  #   Performs column-wise multiplication of selected input columns only.
  # - If condition and input are vectors:
  #   Sets elements in input which are FALSE to 1, and returns it.
  # 
  # Row-wise multiplication means column 1 times column 2 times... e.g. along
  # rows. 
  #
  # Parameters:
  #   condition: a logical vector or a logical matrix. See above.
  #   input: input matrix with data for which to perform calculation.
  # Return: 
  #   - if condition is always FALSE, returns the scalar 1.0
  #   - otherwise a vector with the same length as there are rows in input
  # 
  # Note: Really should be done through R's method thingie...
  #       Maybe once I've read up on it.

  if(!any(condition)) return(1)

  # condition is a matrix
  if(is.matrix(condition)) {  
    # Both are matrices.
    if(!is.matrix(input)) {
      if(ncol(condition) != length(input)) 
        stop("condition does not have as many columns as input has elements.") 
      # Builds input matrix
      input = matrix(input, ncol=length(input), nrow=nrow(condition), byrow=T)
    } else if(any(dim(condition) != dim(input)))
      stop("condition and input have different dimensions.") 

    input[!condition] = 1
    return(prod.matrix(input))
  }

  # condition is a vector and input is a matrix.
  if(is.matrix(input)) {
    if(ncol(input) != length(condition))
      stop("condition does not have as many elements as input has columns.") 

    nbAlleles = sum(condition)
    if(sum(condition) != 1) {
      input = input[, condition, drop=FALSE] 
      input[is.na(input)] = 1 
      return(prod.matrix(input))
    }

    output = input[, condition]
    output[is.na(output)] = 1
    return(output)
  }

  # condition and input are vectors
  if(length(input) != length(condition))
    stop("condition does not have as many elements as input.") 
  output = input
  output[!condition] = 1.0
  output[is.na(output)] = 1.0
  return(output)
}


# Defines prod.matrix from R.
prod.matrix.col.R <- function(x) {
  # Fast column-wise product.
  #
  # Sadly, this is faster than apply(x, 1, prod)  
  y=x[1, ]
  for(i in 2:nrow(x))
  y=y*x[i, ]
  return(y)
}
# Tries and defines a fast prod.matrix if possible.
if(require("inline")) {
  code = "SEXP res;
          int const nrow = INTEGER(GET_DIM(x))[0];
          int const ncol = INTEGER(GET_DIM(x))[1];
          PROTECT(res = allocVector(REALSXP, ncol));
          double *inptr = REAL(x);
          double *const outptr_first = REAL(res);
          double *outptr = outptr_first;
          // Perform multiplication with other columns
          for(int j=1; j < ncol; ++j, ++outptr) {
            *outptr = *inptr; ++inptr;
            if(*outptr == NA_REAL) *outptr = 1.0;
            for(int i=1; i < nrow; ++i, ++inptr)
              if(*inptr != NA_REAL) *outptr *= *inptr;
          }
          UNPROTECT(1);
          return res;"
  prod.matrix.col.C = cfunction(signature(x="array"), code)
  prod.matrix.col = prod.matrix.col.C
} else { prod.matrix.col = prod.matrix.col.R }
prod.matrix.col = prod.matrix.col.R

selective.col.prod <- function(condition, input) {
  # Row-wise product over selected columns or elements.
  #
  # If all(condition == FALSE), returns 1. Otherwise does one of the following.
  #
  # - If condition and input matrices: 
  #   Performs column-wise multiplication of input elements, ignoring those
  #   which are *not* selected by the condition matrix.
  # - If condition is a matrix and input a vector: 
  #   Equivalent to creating a matrix whith the same dimensions as condition
  #   from the input, where row ``i`` is ``input[condition[i, ]]``, then apply
  #   product.
  # - If condition is a vector and input is a matrix: 
  #   Performs column-wise multiplication of selected input columns only.
  # - If condition and input are vectors:
  #   Sets elements in input which are FALSE to 1, and returns it.
  # 
  # Row-wise multiplication means column 1 times column 2 times... e.g. along
  # rows. 
  #
  # Parameters:
  #   condition: a logical vector or a logical matrix. See above.
  #   input: input matrix with data for which to perform calculation.
  # Return: 
  #   - if condition is always FALSE, returns the scalar 1.0
  #   - otherwise a vector with the same length as there are rows in input
  # 
  # Note: Really should be done through R's method thingie...
  #       Maybe once I've read up on it.

  if(!any(condition)) return(1)

  # condition is a matrix
  if(is.matrix(condition)) {  
    # but input is a vector. Then build a matrix from condition and input.
    if(!is.matrix(input)) {
      if(nrow(condition) != length(input)) 
        stop("condition does not have as many columns as input has elements.") 
      # Builds input matrix
      input = matrix(input, nrow=length(input), ncol=ncol(condition))
    } else if(any(dim(condition) != dim(input)))
      stop("condition and input have different dimensions.") 

    input[!condition] = 1
    return(prod.matrix.col(input))
  }

  # condition is a vector and input is a matrix.
  if(is.matrix(input)) {
    if(nrow(input) != length(condition))
      stop("condition does not have as many elements as input has columns.") 

    nbAlleles = sum(condition)
    if(nbAlleles != 1) {
      input = input[condition, , drop=FALSE] 
      input[is.na(input)] = 1 
      return(prod.matrix.col(input))
    }

    output = input[condition, ]
    output[is.na(output)] = 1
    return(output)
  }

  # condition and input are vectors
  if(length(input) != length(condition))
    stop("condition does not have as many elements as input.") 
  output = input
  output[!condition] = 1.0
  output[is.na(output)] = 1.0
  return(output)
}


