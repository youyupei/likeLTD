dropin.probabilities.R <- function(res, freqMat, csp, unc, zeroAll, rate) {
  # Per locus/per replica objective function.
  #
  # Mostly, this separation makes it easier to benchmark different
  # implementation for speed. This particular implementation is the R fallback
  # if inline package is not present.
  selective.col.prod(csp & zeroAll, rate * freqMat)             *
  selective.col.prod(!csp & !unc & zeroAll, 1 - rate * freqMat) *
  res
}
dropin.probabilities  = dropin.probabilities.R
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
  .dropin.probabilities.C = cfunction(sig, code)
  dropin.probabilities.C <- function(res, freqMat, csp, unc, zeroAll, rate) {
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
    .dropin.probabilities.C(res, freqMat, t(!csp & !unc & zeroAll),
                            t(csp != 0 & zeroAll), rate)
  }
  dropin.probabilities <- dropin.probabilities.C
}

dropout.probabilities.R <- function(res, vDoseDropout, csp, unc, ...) {
  # Per locus/per replica objective function.
  #
  # Mostly, this separation makes it easier to benchmark different
  # implementation for speed. This particular implementation is the R fallback
  # if inline package is not present.
  selective.col.prod(!csp & !unc, vDoseDropout)   *
  selective.col.prod(csp != 0, 1 - vDoseDropout)  *
  res
}
dropout.probabilities = dropout.probabilities.R
if(require("inline")) {
  code = "int const nrow = INTEGER(GET_DIM(vDoseDropout))[0];
          int const ncol = INTEGER(GET_DIM(vDoseDropout))[1];
          int * const condA_first = LOGICAL(condA);
          int * const condB_first = LOGICAL(condB);
          int * const condA_end = condA_first + nrow;
          int *zero_ptr = LOGICAL(zeroAll);
          double *vdose_ptr = REAL(vDoseDropout);
          double *out_ptr = REAL(input);
          // Perform same operation for each column
          for(int j=0; j < ncol; ++j, ++out_ptr) 
          {
            int * condA_ptr = condA_first;
            int * condB_ptr = condB_first;
            for(; condA_ptr != condA_end;
                ++vdose_ptr, ++condA_ptr, ++condB_ptr, ++zero_ptr)
            {
              if(*zero_ptr) continue;
              if(*condA_ptr) *out_ptr *= *vdose_ptr;
              else if(*condB_ptr) *out_ptr *= 1.0 - *vdose_ptr;
            }
          }
          return R_NilValue;"
  sig = signature(input="vector", vDoseDropout="array", condA="vector",
                  condB="vector", zeroAll="array")
  .dropout.probabilities.C = cfunction(sig, code)
  dropout.probabilities.C <- function(res, vDoseDropout, csp, unc, zeroAll) {
    if(length(res) != ncol(vDoseDropout))
      stop("output vector and vDoseDropout have incompatible sizes.")
    if(length(csp) != nrow(vDoseDropout))
      stop("csp and vDoseDropout have incompatible sizes.")
    if(length(unc) != nrow(vDoseDropout))
      stop("unc and vDoseDropout have incompatible sizes.")
    if(any(dim(vDoseDropout) != dim(zeroAll)))
      stop("expected vDoseDropout and zeroAll to match.")
    .dropout.probabilities.C(res, vDoseDropout, !csp & !unc, csp != 0, zeroAll)
    res
  }
  dropout.probabilities <- dropout.probabilities.C
}
dropout.probabilities <- dropout.probabilities.C
dropin.probabilities <- dropin.probabilities.R
