impl_objective.dropin.R <- function(res, freqMat, csp, unc, zeroAll, rate) {
  # Per locus/per replica objective function.
  #
  # Mostly, this separation makes it easier to benchmark different
  # implementation for speed. This particular implementation is the R fallback
  # if inline package is not present.
  selective.row.prod(t(csp & zeroAll), rate * freqMat)             *
  selective.row.prod(t(!csp & !unc & zeroAll), 1 - rate * freqMat) *
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
  impl_objective.dropin.C <- function(res, freqMat, csp, unc, zeroAll, rate) {
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
                             t(csp != 0 & zeroAll), rate)
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
if(require("Rcpp")) {
  code = "NumericVector output = NumericVector(input).size();
          NumericMatrix const vDose(vDose);
          LogicalVector const conda(condA);
          LogicalVector const condb(condB);
          for(int i=0; i < output.size(); ++i)
            for(int j=0; j < conda.size(); ++j)
              if(not std::isnan(vDose(i, j))) {
                if(conda(j)) output(i, j) *= vDose(i, j);
                else if(condb(j)) output(i, j) *= 1.0 - vDose(i, j);
              }
          return output;"
  sig = signature(input="vector", vDoseDropout="array", condA="vector",
                  condB="vector")   
  .impl_objective.dropout.cxx = cxxfunction(sig, code,
                                            includes="#include<cmath>",
                                            plugin="Rcpp", )
  impl_objective.dropout.cxx <- function(res, vDoseDropout, csp, unc) {
    if(length(res) != nrow(vDoseDropout))
      stop("output vector and vDoseDropout have incompatible sizes.")
    if(length(csp) != ncol(vDoseDropout))
      stop("csp and vDoseDropout have incompatible sizes.")
    if(length(unc) != ncol(vDoseDropout))
      stop("unc and vDoseDropout have incompatible sizes.")
    .impl_objective.dropout.cxx(res, vDoseDropout, !csp & !unc, csp != 0)
  }
  impl_objective.dropout <- impl_objective.dropout.cxx
}
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
impl_objective.dropout <- impl_objective.dropout.R
impl_objective.dropin <- impl_objective.dropin.R
