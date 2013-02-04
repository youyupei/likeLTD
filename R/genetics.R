
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
  queried = known.alleles(name, refData)
  # Checks if there are any dropout for specific locus
  for(locus in names(cprofs)) 
    for(replicate in cprofs[[locus]]) {
      # Replicate could be empty or invalid in some way.
      # In that case, we shouldn't call missing alleles as missing. 
      validRep = !any(sapply(replicate$csp, is.na))
      if(validRep && any(!queried[[locus]] %in% replicate$csp)) return(TRUE)
    }
  return(FALSE)
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
  presence.per.locus <- function(frq.locus, cprofs.locus, ..., type="csp") {
    # Same as presence.matrices but for a single locus.
    
    # Check if valid locus first.
    invalidCSP = unlist(sapply(cprofs.locus, function(n) is.na(n$csp)))
    if( any(invalidCSP) ) return(NULL)

    # count: creates a row of the matrix.
    #        A row consists of 0 and 1, where 0 indicates that the allele is
    #        not present in the CSP. It is also possible to add other list of
    #        alleles via the ellipsis. In that case, a row consists of integers
    #        indicating the number of lists in which the allele is found. 
    rep = row.names(frq.locus)
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
    result <- sapply(cprofs.locus, count) 
    return(t(matrix(result, ncol=length(cprofs.locus))))
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

  adjust.per.locus <- function(frq.locus, q.loc, adj=1, fst=0.02) {
    # Applies operations to a single locus.
    homozygote <- as.integer(q.loc[1] == q.loc[2])
    frq.locus[q.loc, 1] <- frq.locus[q.loc, 1] + adj * (1 + homozygote)
    #  shared ancestry adjustements.
    frq.locus[, 1] <- frq.locus[, 1] / sum(frq.locus[, 1]) * (1 - fst) / (1 + fst) 
    frq.locus[q.loc, 1] <- frq.locus[q.loc, 1] + fst / (1 + fst) * (1 + homozygote)
    return(frq.locus)
  }
  return(mapply(adjust.per.locus, frequencies, queriedAlleles, adj=adj,
                fst=fst))
}
