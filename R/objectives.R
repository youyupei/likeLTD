empty.alleles = function(genotypes, dropoutPresence, nUnknowns) {
  # Mask over genotype matrix indicating empty/null elements.
  #
  # Parameters:
  #   genotypes: Matrix of indices into dosage matrix.
  #   dropoutPresence: Profiles with dropout.
  # Returns: A matrix indicating which alleles will be NA in dosage.
  if(!is.matrix(dropoutPresence)) stop("Expected matrix in input.")
  if(ncol(genotypes) == 0) return(matrix(nrow=nrow(dropoutPresence), ncol=0))
  knownZero = rowSums(dropoutPresence) == 0
  if(nUnknowns == 0) 
    return(matrix(knownZero, nrow=length(knownZero), ncol=ncol(genotypes)))
  iota = 1:length(knownZero)
  apply(genotypes, 2, function(n) (!iota %in% n) & knownZero) 
}

relatedness.factors <- function(genotypes, alleleDb, queriedProfile,
                                relatedness=c(0, 0)) {
  # Creates relatedness factor 
  # 
  # Parameters:
  #  genotypes: 
  if(all(abs(relatedness) < 1e-8)) return(1)
  indices = c(which(rownames(alleleDb) %in% queriedProfile[[1]]), 
              which(rownames(alleleDb) %in% queriedProfile[[2]]) )
  defactored = function(n) {
    result = 1.0 - sum(relatedness)
    hasFirst = n %in% indices[[1]]
    hasSecond = n %in% indices[[2]]
    if(any(hasFirst)) {
      factor = relatedness[[1]] * 0.5 / alleleDb[indices[[1]], 1]
      if(!all(hasFirst)) factor = factor * 0.5
      result = result + factor
    }
    if(any(hasSecond)) {
      factor = relatedness[[1]] * 0.5 / alleleDb[indices[[2]], 1]
      if(!all(hasSecond)) factor = factor * 0.5
      result = result + factor
    }
    if(any(hasFirst) & any(hasSecond)) {
      factor = relatedness[[2]] / alleleDb[indices[[1]], 1] /
               alleleDb[indices[[2]], 1] 
      if(indices[[1]] != indices[[2]]) factor = factor * 0.5
      result = result + factor
    }
    return(result)
  }
  apply(genotypes[1:2, ], 2, defactored) 
}

genotype.factors <- function(genotypes, alleleDb, nUnknowns, doDropin,
                             queriedProfile, relatedness=c(0, 0)) {
  # Fraction cum heterozygote factors.
  #

  # First add standard bits.
  if(nUnknowns == 0 & !doDropin) het = matrix(1, ncol=1, nrow=nrow(alleleDb))
  else {
    n = nrow(genotypes) / 2
    odd = 2*1:n - 1; even = 2*1:n
    het <- 1 + matrix(genotypes[odd, ] < genotypes[even, ], nrow=length(odd))
    het <- apply(het, 2, prod)
  }
  fractions <- matrix(alleleDb[genotypes, 1], ncol=ncol(genotypes))
  result <- apply(fractions, 2, prod) * het
  result * relatedness.factors(genotypes, alleleDb, queriedProfile,
                               relatedness)
}

likelihood.constructs.per.locus = function(scenario) {
  # Creates the locus-specific data needed by each likehood function.
  #
  # Parameters:
  #   scenario: A scenario, for instance one returned by
  #             prosecution.scenario(...) or defense.scenario(...)
  alleles = rownames(scenario$alleleDb)
  if(is.null(alleles)) stop("Could not figure out alleles names.")
  alleles.vector = function(n) alleles %in% unlist(n)
  cspPresence     = apply(scenario$cspProfile, 1, alleles.vector)
  dropoutPresence = apply(scenario$dropoutProfs, 1, alleles.vector)
  uncPresence     = apply(scenario$uncProf, 1, alleles.vector)
  if(!is.matrix(cspPresence))
    cspPresence = matrix(nrow=0, ncol=length(alleles))
  if(!is.matrix(dropoutPresence))
    dropoutPresence = matrix(nrow=0, ncol=length(alleles))
  if(!is.matrix(uncPresence))
    uncPresence = matrix(nrow=0, ncol=length(alleles))

  missingReps = apply(scenario$cspProfile, 1, is.na)

  genotypes <- possible.genotypes(cspPresence, dropoutPresence, missingReps,
                                  alleles, scenario$nUnknowns,
                                  scenario$doDropin)
  zeroAll = empty.alleles(genotypes, dropoutPresence, scenario$nUnknowns) 
  factors = genotype.factors(genotypes, scenario$alleleDb, scenario$nUnknowns,
                             scenario$doDropin, scenario$queriedProfile,
                             scenario$relatedness) 
  freqMat = matrix(scenario$alleleDb[, 1], ncol=ncol(genotypes),
                   nrow=nrow(scenario$alleleDb))
  
  fragLengths = matrix(scenario$alleleDb[genotypes, 2], ncol=ncol(genotypes))

  list(cspPresence=cspPresence, dropoutPresence=dropoutPresence,
       uncPresence=uncPresence, missingReps=missingReps,
       genotypes=genotypes, zeroAll=zeroAll, factors=factors,
       freqMat=freqMat, fragLengths=fragLengths)
}

create.likelihood.per.locus <- function(scenario, addAttr=FALSE) {
  # Creates a likelyhood function for a given scenario and locus
  #
  # A scenario is given by the number of unknown contributors, whether to model
  # dropin, so on and so forth.

  cons = likelihood.constructs.per.locus(scenario)

  result.function <- function(rcont, degradation, localAdjustment,
                              tvedebrink, dropout, dropin, ...) {
    # Likelyhood function for a given scenario and locus
    #
    # This function is specific to the scenario for which it was created.
    # It hides everything except the nuisance parameters over which to
    # optimize.
    #
    # Parameters:
    #   rcont: relative contribution from each profiled individual and unknown
    #          contributor in this scenario.
    #   degradation: relative degradation from each profiled individual in this
    #                scenario.
    #   localAdjustment: a scalar floating point value.
    #   tvedebrink: a scalar floating point value.
    #   dropout: the dropout rate for each replicate.
    #   ...: Any other parameter, e.g. for penalty functions. 
    #        These parameters are ignored here.
    # Returns: A scalar value giving the likelihood for this locus and
    #          scenario.
    if(length(rcont) != scenario$nUnknowns + ncol(cons$dropoutPresence))
      stop(sprintf("rcont should be %d long.",
                   scenario$nUnknowns + ncol(cons$dropoutPresence)))
    if(length(degradation) != scenario$nUnknowns + ncol(cons$dropoutPresence))
      stop(sprintf("degradation should be %d long.",
                   scenario$nUnknowns + ncol(cons$dropoutPresence)))
    if(length(dropout) != length(cons$missingReps))
      stop(sprintf("dropout should be %d long.", ncol(cons$dropoutPresence)))

    allEPG <- all.epg.per.locus(rcont, degradation, cons$dropoutPresence,
                                scenario$alleleDb[, 2], cons$fragLengths,
                                cons$genotypes, scenario$nUnknowns > 0)
    allEPG = (allEPG * localAdjustment)^tvedebrink

    # res: (Partial) Likelihood per allele.
    res = array(1, length(cons$factors))
    # Loop over replicates.
    for(i in 1:length(cons$missingReps)) {
      if(!cons$missingReps[i]) {
        csp = cons$cspPresence[, i]
        unc = cons$uncPresence[, i]

        vDoseDropout = allEPG * dropout[i]
        vDoseDropout = vDoseDropout / (vDoseDropout + 1 - dropout[i])

        res = dropout.probabilities(res, vDoseDropout, csp, unc, cons$zeroAll)
        if(scenario$doDropin) 
          res = dropin.probabilities(res, cons$freqMat, csp, unc, cons$zeroAll,
                                     dropin * (1 - dropout[i]) )
      } # End of if(any(csp)) 
    } # End of loop over replicates.

    # Figure out likelihood for good and return.
    return(sum(res * cons$factors))
  }

  if(addAttr) {
    attr(result.function, "scenario") <- scenario
    attr(result.function, "constructs") <- cons
  }
  result.function
}


penalties <- function(rcont, degradation, localAdjustment, tvedebrink, dropout,
                      dropin, localAdjPenalty=50, dropinPenalty=2,
                      degradationPenalty=50, bemn=-4.35, besd=0.38, ...)
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
  # Normalizes by number of loci so product of penalties same as in old code.
  # Some penalties are per locus and other for all locus. Hence this bit of run
  # around.
  normalization = 1.0 / max(length(localAdjustment), 1)

  if(!missing(rcont) & !is.null(rcont))
    result = result * exp(-dropin * dropinPenalty * normalization)
  if(!missing(degradation) & !is.null(degradation))
    result = result * exp(-sum(degradation) * degradationPenalty 
                           * normalization )
  if(!missing(tvedebrink) & !is.null(tvedebrink))
    result = result * dnorm(tvedebrink, bemn, besd) ^ normalization
  if(!missing(localAdjustment) & !is.null(localAdjustment))
    result = result * dgamma(localAdjustment, localAdjPenalty, localAdjPenalty)

  return(result)
}

create.likelihood.vectors <- function(scenario, addAttr=FALSE) {
  # Creates a likelihood function from the input scenario.
  #
  # Likelihood function returns a list with two objects: an array of likelihood
  # per locus, an array of penalties per locus.
  #
  # .. seealso:: create.likelihood, create.likelihood.log
  #
  # Parameters
  #   scenario: A list containing all the items necessary to building the
  #             objective function. It is likely the result of a call to
  #             defense.scenario or prosecution.scenario.
  #   addAttr: If TRUE, adds some attribute to the objective function.
  #            This is mostly for debug purposes. It makes it easier to known
  #            what parameters were use to build the objective function.

  locusCentric = transform.to.locus.centric(scenario)
  functions <- mapply(create.likelihood.per.locus, locusCentric,
                      MoreArgs=list(addAttr=addAttr))

  likelihood.vectors <- function(rcont, degradation, localAdjustment,
                                 tvedebrink, dropout, dropin,
                                 localAdjPenalty=50, dropinPenalty=2,
                                 degradationPenalty=50, bemn=-4.35,
                                 besd=0.38, ...) {
    # Call each and every function in the array.
    arguments = list(rcont=rcont, degradation=degradation,
                     tvedebrink=tvedebrink, dropout=dropout,
                     dropinPenalty=dropinPenalty,
                     degradationPenalty=degradationPenalty, bemn=bemn,
                     besd=besd, dropin=dropin)
    callme <- function(objective, adj) {
      args = append(arguments, list(localAdjustment=adj))
      do.call(objective, args)
    }
    if(length(localAdjustment) == 1)
      localAdjustment = rep(localAdjustment, length(functions))
    objectives = mapply(callme, functions, localAdjustment)
    arguments = append(arguments, list(localAdjustment=localAdjustment))
    arguments = append(arguments, list(...))
    pens <- do.call(penalties, arguments)
    list(objectives=objectives, penalties=pens)
  }
  if(addAttr) {
    attr(likelihood.vectors, "scenario") <- scenario
    attr(likelihood.vectors, "functions") <- functions
  }
  return(likelihood.vectors)
}


create.likelihood <- function(scenario, addAttr=FALSE) {
  # Creates a likelihood function from the input scenario.
  #
  # .. seealso:: create.likelihood.vectors, create.likelihood.log

  vecfunc <- create.likelihood.vectors(scenario, addAttr)

  likelihood.scalar <- function(...) { 
    result <- vecfunc(...) 
    prod(result$objectives * result$penalties)
  }

  attributes(likelihood.scalar) <- attributes(vecfunc)
  return(likelihood.scalar)
}

create.likelihood.log <- function(scenario, addAttr=FALSE) {
  # Creates a likelihood function from the input scenario.
  #
  # .. seealso::
  #   
  #   create.likelihood, create.likelihood.vectors
  vecfunc <- create.likelihood.vectors(scenario, addAttr)
  likelihood.log <- function(...) { 
    result <- vecfunc(...) 
    sum(log(result$objectives)) # + sum(log(result$penalties))
  }
  attributes(likelihood.log) <- attributes(vecfunc)
  return(likelihood.log)
}
