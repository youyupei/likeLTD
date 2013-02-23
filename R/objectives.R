impl_objective.dropin.R <- function(res, freqMat, csp, unc, zeroAll, dropout) {
  # Per locus/per replica objective function.
  #
  # Mostly, this separation makes it easier to benchmark different
  # implementation for speed. This particular implementation is the R fallback
  # if inline package is not present.
  selective.row.prod(t(csp & zeroAll), (1 - dropout) * freqMat)             *
  selective.row.prod(t(!csp & !unc & zeroAll), 1 - (1 - dropout) * freqMat) *
  res
}
impl_objective.dropin  = impl_objective.dropin.R
if(require("inline")) {
  code = "SEXP result;
          PROTECT(result = Rf_duplicate(input));
          int const nrow = INTEGER(GET_DIM(condA))[0];
          int const ncol = INTEGER(GET_DIM(condA))[1];
          int *condA_ptr = LOGICAL(condA);
          int *condB_ptr = LOGICAL(condB);
          double *freqmat_ptr = REAL(freqMat);
          double *const outptr_first = REAL(result);
          double const rate = 1 - *REAL(dropinRate);
          for(int j=0; j < ncol; ++j) 
          {
            double * outptr = outptr_first;
            for(int i=0; i < nrow; ++i, ++condA_ptr, ++condB_ptr,
                                   ++freqmat_ptr, ++outptr) 
              if(*condA_ptr) *outptr *= 1.0 - rate * (*freqmat_ptr);
              else if(*condB_ptr) *outptr *= rate * (*freqmat_ptr);
          }
          UNPROTECT(1);
          return result;"
  sig = signature(input="vector", freqMat="array", condA="vector",
                  condB="vector", dropinRate="real")
  .impl_objective.dropin.C = cfunction(sig, code)
  impl_objective.dropin.C <- function(res, freqMat, csp, unc, zeroAll,
                                      dropout) {
    if(length(res) != nrow(freqMat))
      stop("output vector and freqMat have incompatible sizes.")
    if(length(csp) != ncol(freqMat))
      stop("csp and vDoseDropout have incompatible sizes.")
    if(length(unc) != length(csp))
      stop("unc and csp have incompatible sizes.")
    if(length(unc) != nrow(zeroAll))
      stop("unc and zeroAll have incompatible sizes.")
    if(nrow(freqMat) != ncol(zeroAll))
      stop("unc and zeroAll have incompatible sizes.")
    .impl_objective.dropin.C(res, freqMat, t(!csp & !unc & zeroAll),
                             t(csp != 0 & zeroAll), dropout)
  }
  impl_objective.dropin <- impl_objective.dropin.C
}

impl_objective.dropout.R <- function(res, vDoseDropout, csp, unc) {
  # Per locus/per replica objective function.
  #
  # Mostly, this separation makes it easier to benchmark different
  # implementation for speed. This particular implementation is the R fallback
  # if inline package is not present.
  selective.row.prod(!csp & !unc, vDoseDropout)   *
  selective.row.prod(csp != 0, 1 - vDoseDropout)  *
  res
}
impl_objective.dropout = impl_objective.dropout.R
if(require("inline")) {
  code = "SEXP result;
          PROTECT(result = Rf_duplicate(input));
          int const nrow = INTEGER(GET_DIM(vDoseDropout))[0];
          int const ncol = INTEGER(GET_DIM(vDoseDropout))[1];
          int *condA_ptr = LOGICAL(condA);
          int *condB_ptr = LOGICAL(condB);
          double * vdose_ptr = REAL(vDoseDropout);
          double *const outptr_first = REAL(result);
          // Perform same operation for each column
           for(int j=0; j < ncol; ++j, ++condA_ptr, ++condB_ptr) 
             if(*condA_ptr or condB_ptr) {
               double *outptr = outptr_first;
               if(*condA_ptr) {
                 for(int i=0; i < nrow; ++i, ++vdose_ptr, ++outptr)
                   if(not isnan(*vdose_ptr)) *outptr *= *vdose_ptr;
               } else { // If we are here, condB_ptr is true.
                 for(int i=0; i < nrow; ++i, ++vdose_ptr, ++outptr)
                   if(not isnan(*vdose_ptr)) *outptr *= 1.0 - *vdose_ptr;
               } 
             }
          UNPROTECT(1);
          return result;"
  sig = signature(input="vector", vDoseDropout="array", condA="vector",
                  condB="vector")
  .impl_objective.dropout.C = cfunction(sig, code)
  impl_objective.dropout.C <- function(res, vDoseDropout, csp, unc) {
    if(length(res) != nrow(vDoseDropout))
      stop("output vector and vDoseDropout have incompatible sizes.")
    if(length(csp) != ncol(vDoseDropout))
      stop("csp and vDoseDropout have incompatible sizes.")
    if(length(unc) != ncol(vDoseDropout))
      stop("unc and vDoseDropout have incompatible sizes.")
    .impl_objective.dropout.C(res, vDoseDropout, !csp & !unc, csp != 0)
  }
  impl_objective.dropout <- impl_objective.dropout.C
}

create.likelihood.per.locus <- function(queriedPresence, profPresence,
                                        cspPresence, uncPresence, missingReps,
                                        alleleDb, nUnknowns, doDropin) {
  # Creates a likelyhood function for a given scenario and locus
  #
  # A scenario is given by the number of unknown contributors, whether to model
  # dropin, so on and so forth.

  #############################################################################
  ##############################  PREPARATION  ################################
  #############################################################################
  # All possible sets of alleles from unknown contributors within this
  # scenario.
  allProfiles <- possible.profiles(queriedPresence, cspPresence, profPresence,
                                   missingReps, row.names(alleleDb),
                                   nUnknowns, doDropin)
  
  # FragLengths: lengths of each short tandem repeat for each profile.
  FragLengths = matrix(alleleDb[allProfiles, 2],
                       ncol=nrow(allProfiles), byrow=T)

  # Mask for empty alleles in each epg 
  knownZero = !( colSums(profPresence) > 0 )
  if(nUnknowns == 0)
    zeroAll = matrix(knownZero, ncol=nrow(allProfiles), nrow=nrow(alleleDb))
  else
    zeroAll = apply(allProfiles, 1, 
                    function(n) { 
                      result <- knownZero; result[n] = FALSE; result } )

  # Frequencies as a matrix, for fater computing
  freqMat = matrix(alleleDb[, 1], nrow=nrow(allProfiles),
                   ncol=length(alleleDb[, 1]), byrow=T)
  # nrep: Number of replicates. 
  nrep = nrow(cspPresence)

  # Prepare heterozygote adjustment.
  if(nUnknowns > 0) {
    het <- 1 + (allProfiles[, 2*(1:nUnknowns)-1] <
                  allProfiles[, 2*(1:nUnknowns)])
  } else het <- 1 + (allProfiles[, 1] < allProfiles[, 2])
  if (nUnknowns > 1) het <- apply(het, 1, prod, na.rm=T)

  # Prepare allele fractions 
  fraction <- matrix(alleleDb[allProfiles, 1], ncol=ncol(allProfiles))
  fraction <- apply(fraction, 1, prod, na.rm=T)

  #############################################################################
  ####################  PER LOCUS OBJECTIVE FUNCTION  #########################
  #############################################################################
  result.function <- function(rcont, degradation, localAdjustment,
                              tvedebrink, dropout, ...) {
    # Likelyhood function for a given scenario and locus
    #
    # This function is specific to the scenario for which it was created.
    # It hides everything except the nuisance parameters over which to
    # optimize.
    #
    # Parameters:
    #   rcont: relative contribution from each profiled individual in this
    #          scenario.
    #   degradation: relative degradation from each profiled individual in this
    #          scenario.
    #   localAdjustment: a scalar floating point value.
    #   tvedebrink: a scalar floating point value.
    #   dropout: the dropout rate for each replicate.
    #   ...: Any other parameter, e.g. for penalty functions. 
    #        These parameters are ignored.
    # Returns: A scalar value giving the likelihood for this locus and
    #          scenario.
    allEPG <- all.epg.per.locus(rcont, degradation, profPresence,
                                alleleDb[, 2], FragLengths, allProfiles,
                                nUnknowns > 0)
    allEPG = t(allEPG * localAdjustment)^tvedebrink


    # res: (Partial) Likelihood per allele.
    res = array(1, length(fraction))
    # Loop over replicates.
    for(i in 1: nrep) {
      if(!missingReps[i]) {
        csp = cspPresence[i, ]
        unc = uncPresence[i, ]

        vDoseDropout = allEPG * dropout[i]
        vDoseDropout = vDoseDropout / (vDoseDropout + 1 - dropout[i])

        res = impl_objective.dropout(res, vDoseDropout, csp, unc)
        if(doDropin) 
          res = impl_objective.dropin(res, freqMat, csp, unc, zeroAll,
                                      dropout[i])
      } # End of if(any(csp)) 
    } # End of loop over replicates.

    # Figure out likelihood for good and return.
    return(sum(res * fraction * het))
  }

  ############################################################################# 
  ####################  RETURNING OBJECTIVE FUNCTION  #########################
  ############################################################################# 
  attr(result.function, "nUnknowns")       <- nUnknowns
  attr(result.function, "doDropin")        <- doDropin
  attr(result.function, "alleleDb")        <- alleleDb
  attr(result.function, "queriedPresence") <- queriedPresence
  attr(result.function, "profPresence")    <- profPresence
  attr(result.function, "cspPresence")     <- cspPresence
  attr(result.function, "uncPresence")     <- uncPresence
  attr(result.function, "missingReps")     <- missingReps
  return(result.function)
}


create.likelihood.functions <- function(admin, nUnknowns=0, doDropin=FALSE,
                                        fst=0.02, adj=1.0, ethnic='EA1') {
  # Creates vector of likelihood functions, one per locus.
  #
  # This function returns an array of functions. The point is to make it easier
  # to call each locus separately, if necessary. 
  # 
  # .. seealso::
  #   
  #   create.likelihood, create.log.likelihood, create.likelihood.vectors

  # For each locus, we need to create the input to create.likelihood.per.locus.
  # This is most likely not the shortest path. It is, however, that one that
  # was implemented first...
  genetics <- pack.genetics.input(admin, unknowns=nUnknowns, dropin=doDropin)
  knownAlleles <- known.alleles( c(genetics$nameQ, genetics$nameK),
                                 genetics$refData )
  alleleDb <- ethnic.database(ethnic, names(genetics$cprofs), genetics$afreq)
  alleleDb <- add.missing.alleles( alleleDb, genetics$cprofs,
                                   knownAlleles[1:2, ] )
  alleleDb <- adjust.frequencies(alleleDb, knownAlleles[1:2, ], adj=adj,
                                 fst=fst)
  knownAllelesNoDropouts = known.without.dropouts( genetics$nameK,
                                                   genetics$refData,
                                                   genetics$cprofs )
  cspPresence = presence.matrices(alleleDb, genetics$cprofs)
  if(length(knownAllelesNoDropouts) == 0) {
    uncPresence = presence.matrices(alleleDb, genetics$cprofs, type="unc")
  } else {
    uncPresence = presence.matrices(alleleDb, genetics$cprofs, knownAllelesNoDropouts,
                                    type="unc")
  }
  presence.per.locus <- function(prof, alleleDb) {
    matrix( sapply(prof, function(n) as.integer(rownames(alleleDb) %in% n)), 
            nrow = length(prof), byrow=T )
  }
  indices = (length(genetics$nameQ) * 2 + 1):nrow(knownAlleles)
  profPresence <- mapply(presence.per.locus,
                         knownAlleles[indices, names(alleleDb)], alleleDb)
  indices = 1:(length(genetics$nameQ) * 2)
  queriedPresence <- mapply(presence.per.locus,
                            knownAlleles[indices, names(alleleDb)], alleleDb)
  missingReps = lapply(genetics$cprofs, missing.csp.per.locus)

  return( mapply(create.likelihood.per.locus, 
                 queriedPresence,
                 profPresence,
                 cspPresence,
                 uncPresence, 
                 missingReps,
                 alleleDb, 
                 MoreArgs = list(nUnknowns=nUnknowns, doDropin=doDropin)) )
}

penalty.localAdj <- function(value, penalty=50) dgamma(value, penalty, penalty)
penalty.dropin      <- function(value, penalty=2)  exp(-value[length(value)] * penalty)
penalty.degradation <- function(value, penalty=50) exp(-sum(value) * penalty)
penalties <- function(localAdjustment, localAdjPenalty=50, rcont,
                      dropinPenalty=2, degradation, degradationPenalty=50,
                      beta, bemn=-4.35, besd=0.38, ...)
{
  # Penalties to apply to the likelihood.
  #
  # The return can be a scalar (no localAdj argument) or a vector (localAdj
  # argument). They can be applied to the likelihood as:
  #
  #   >>> prod(likelihood * penalties(...))
  #
  # Where likelihood is a vector where each element is the likelihood for a
  # given locus.
  #
  # It should be applied to the negative log likelihood as:
  #
  #  >>> sum( -log(likelihood) - log(penalties(...)))
  result = 1
  if(!missing(rcont) & !is.null(rcont))
    result = result * penalty.dropin(rcont, dropinPenalty)
  if(!missing(degradation) & !is.null(degradation))
    result = result * penalty.degradation(degradation, degradationPenalty)
  if(!missing(beta) & !is.null(beta))
    result = result * dnorm(beta, bemn, besd)
  if(!missing(localAdjustment) & !is.null(localAdjustment))
    result = result * penalty.localAdj(localAdjustment, localAdjPenalty)

  return(result)
}

create.likelihood.vectors <- function(admin, nUnknowns=0, doDropin=FALSE,
                                      fst=0.02, adj=1.0, ethnic='EA1') {
  # Creates a likelihood function from the input scenario.
  #
  # Likelihood function returns a list with two objects: an array of likelihood
  # per locus, an array of penalties per locus.
  # .. seealso::
  #   
  #   create.likelihood, create.log.likelihood, create.likelihood.functions

  function.array <- create.likelihood.functions(admin=admin,
                                                nUnknowns=nUnknowns,
                                                doDropin=doDropin, fst=fst,
                                                adj=adj, ethnic=ethnic)

  likelihood.vectors <- function(rcont, degradation, localAdjustment,
                                 tvedebrink, dropout, beta, localAdjPenalty=50,
                                 dropinPenalty=2, degradationPenalty=50,
                                 bemn=-4.35, besd=0.38) {
    # Call each and every function in the array.
    arguments = list(rcont=rcont, degradation=degradation,
                     tvedebrink=tvedebrink, dropout=dropout,
                     beta=beta, dropinPenalty=dropinPenalty,
                     degradationPenalty=degradationPenalty, bemn=bemn,
                     besd=besd)
    callme <- function(objective, adj) {
      args = append(arguments, list(localAdjustment=adj))
      do.call(objective, args)
    }
    if(length(localAdjustment) == 1)
      localAdjustment = rep(localAdjustment, length(function.array))
    objectives = mapply(callme, function.array, localAdjustment)
    pens <- do.call( penalties, 
                     append(arguments, list(localAdjustment=localAdjustment)) )
    list(objectives=objectives, penalties=pens)
  }
  attr(likelihood.vectors, "admin")     <- admin
  attr(likelihood.vectors, "nUnknowns") <- nUnknowns
  attr(likelihood.vectors, "doDropin")  <- doDropin
  attr(likelihood.vectors, "fst")       <- fst
  attr(likelihood.vectors, "adj")       <- adj
  attr(likelihood.vectors, "ethnic")    <- ethnic
  attr(likelihood.vectors, "loci")      <- function.array

  return(likelihood.vectors)
}


create.likelihood <- function(admin, nUnknowns=0, doDropin=FALSE, fst=0.02,
                              adj=1.0, ethnic='EA1') {
  # Creates a likelihood function from the input scenario.
  #
  # .. seealso::
  #   
  #   create.likelihood, create.log.likelihood, create.likelihood.functions

  vecfunc <- create.likelihood.vectors(admin=admin, nUnknowns=nUnknowns,
                                       doDropin=doDropin, fst=fst, adj=adj,
                                       ethnic=ethnic)

  likelihood.scalar <- function(...) { 
    result <- vecfunc(...) 
    prod(result$objectives) # * result$penalties)
  }

  attributes(likelihood.scalar) <- attributes(vecfunc)
  return(likelihood.scalar)
}

create.likelihood.log <- function(admin, nUnknowns=0, doDropin=FALSE, fst=0.02,
                                  adj=1.0, ethnic='EA1') {
  # Creates a likelihood function from the input scenario.
  #
  # .. seealso::
  #   
  #   create.likelihood, create.likelihood.log, create.likelihood.functions
  vecfunc <- create.likelihood.vectors(admin=admin, nUnknowns=nUnknowns,
                                       doDropin=doDropin, fst=fst, adj=adj,
                                       ethnic=ethnic)
  likelihood.log <- function(...) { 
    result <- vecfunc(...) 
    sum(log(result$objectives)) # + sum(log(result$penalties))
  }

  attributes(likelihood.log) <- attributes(vecfunc)
  return(likelihood.log)
}
