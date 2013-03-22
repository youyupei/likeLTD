probabilities.function <- function(hypothesis, cons, doR=FALSE) {
  # Creates a probability function.
  #
  # The probability function computes all probabilities associated with a
  # specific set of genotypes. 
  #
  # In practice, this function makes is easy to substitute C vs R
  # implementations, as well as specialize the C implementations. 

  if((!hypothesis$doDropin) && !doR) {
    probabilities.no.dropin.C <- function(res, vDoseDropout, csp, unc, ...) {
        if(length(res) != ncol(vDoseDropout))
          stop("output vector and vDoseDropout have incompatible sizes.")
        if(length(csp) != nrow(vDoseDropout))
          stop("csp and vDoseDropout have incompatible sizes.")
        if(length(unc) != nrow(vDoseDropout))
          stop("unc and vDoseDropout have incompatible sizes.")
        if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
          stop("expected vDoseDropout and zeroAll to match.")
        .Call(.cpp.probabilitiesNoDropin, res, vDoseDropout, !csp & !unc, csp,
              cons$zeroAll, PACKAGE="likeLTD") 
    }
    return(probabilities.no.dropin.C) 
  } else if(hypothesis$doDropin && !doR) {
    probabilities.with.dropin.C <- function(res, vDoseDropout, csp, unc, rate) {
       if(length(res) != ncol(vDoseDropout))
         stop("output vector and vDoseDropout have incompatible sizes.")
       if(length(csp) != nrow(vDoseDropout))
         stop("csp and vDoseDropout have incompatible sizes.")
       if(length(unc) != nrow(vDoseDropout))
         stop("unc and vDoseDropout have incompatible sizes.")
       if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
         stop("expected vDoseDropout and zeroAll to match.")
       if(any(dim(cons$freqMat) != dim(cons$zeroAll)))
         stop("expected vDoseDropout and zeroAll to match.")
       .Call(.cpp.probabilitiesWithDropin, res, vDoseDropout, !csp & !unc, csp,
             cons$zeroAll, cons$freqMat, rate, PACKAGE="likeLTD")
     }
    return(probabilities.with.dropin.C)
  }
  # Otherwise, return an R function. 
  function(res, vDoseDropout, csp, unc, rate) {
     res <- res                                             *
            selective.col.prod(!csp & !unc, vDoseDropout)   *
            selective.col.prod(csp, 1 - vDoseDropout)  
    if(hypothesis$doDropin) 
      res <- res                                                           *
             selective.col.prod(csp & cons$zeroAll, rate * cons$freqMat)   *
             selective.col.prod(!csp & !unc & cons$zeroAll,
                                1 - rate * cons$freqMat)
    res
  }
}

all.genotypes.per.locus.C <- function(nAlleles, nContrib=1) {
  if(nAlleles < 0) stop("Negative or null number of alleles.")
  if(nContrib < 0) stop("Negative number of contributors.")
  if(nContrib == 0 || nAlleles == 0) return(matrix(nrow=0, ncol=0))
  nContrib = as.integer(nContrib)
  a = t(combinations(nAlleles, 2, repeats.allowed=TRUE))
  if(is.null(nContrib) || is.null(a)) stop("something went very wrong.")
  .Call(.cpp.allGenotypesPerLocus, nContrib, a, PACKAGE="likeLTD")
}
all.genotypes.per.locus = all.genotypes.per.locus.C
