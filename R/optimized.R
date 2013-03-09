# compile.C <- function(sig, code) cfunction(sig, code)
# if(require("inline")) {
#   openmp.params <- function() {
#     if(!require("inline")) return(NULL)
#     # Try fopenmp
#     result = list(includes=c("#include <omp.h>"), 
#                   cxxargs=c("-fopenmp"),
#                   libargs=c("-fopenmp"))
#     code = "#ifdef _OPENMP
#               return 0;
#             #else
#               breaks_on_purpose
#             #endif"
#     error = FALSE
#     error.function = function(e) { error <<- TRUE }
#     tryCatch(cfunction(signature(), code, includes=result$includes,
#                        cxxargs=result$cxxargs, libargs=result$libargs),
#              error=error.function)
#     if(!error) return(result)
#     return(NULL)
#   }
#   openmpParams <- openmp.params()
#   if(!is.null(openmpParams)) {
#     compile.C <- function(sig, code) {
#       return(cfunction(sig, code, 
#                        includes=c(openmpParams$includes,
#                                   "#define CSTACK_DEFNS 7 ",
#                                   "#include <Rinterface.h>"),
#                        cxxargs=openmpParams$cxxargs,
#                        libargs=openmpParams$libargs))
#     }
#   }
# } 
# if(require("inline")) {
#   openmp.diagnostic <- compile.C(
#     signature(),
#     code = "SEXP result;
#             #ifdef _OPENMP
#               int nbthreads;
#               #pragma omp parallel shared(nbthreads)
#               { 
#                 #pragma omp single
#                 nbthreads = omp_get_num_threads();
#               }
#             #else
#               int const nbthreads = 0;
#             #endif
#             PROTECT(result = NEW_INTEGER(1));
#             INTEGER_POINTER(result)[0] = nbthreads;
#             UNPROTECT(1);
#             return result;"
#   )
#   .probabilities.no.dropin.C <- compile.C(
#     signature(input="vector", vDoseDropout="array", condA="vector",
#               condB="vector", zeroAll="array"),
#     code = "#ifdef _OPENMP
#               uintptr_t const oldstack = R_CStackLimit;
#               R_CStackLimit = (uintptr_t) - 1;
#             #endif
#             int const nrow = INTEGER(GET_DIM(vDoseDropout))[0];
#             int const ncol = INTEGER(GET_DIM(vDoseDropout))[1];
#             int const * const condA_first  = LOGICAL(condA);
#             int const * const condB_first  = LOGICAL(condB);
#             int const * const zero_ptr     = LOGICAL(zeroAll);
#             double const * const vdose_ptr = REAL(vDoseDropout);
#             double * const out_ptr         = REAL(input);
#             
#             #pragma omp parallel for schedule(static, 1000)
#             for(int j=0; j < ncol; ++j) 
#             {
#               int const v = j * nrow;
#               int const *condA_ptr = condA_first;
#               int const *condB_ptr = condB_first;
#               for(int i=0; i < nrow; ++i, ++condA_ptr, ++condB_ptr) 
#               {
#                 if(*(zero_ptr + v + i)) continue;
#                 if(*condA_ptr)      
#                   *(out_ptr + j) *= *(vdose_ptr + v + i);
#                 else if(*condB_ptr) 
#                   *(out_ptr + j) *= 1.0 - *(vdose_ptr + v + i);
#               }
#             }
#             #ifdef _OPENMP
#                R_CStackLimit = oldstack;
#             #endif
#             
#             return R_NilValue;"
#   )

#   .probabilities.with.dropin.C <- compile.C(
#     signature(input="vector", vDoseDropout="array", condA="vector",
#               condB="vector", zeroAll="array", freqMat="array", dropin="real"),
#     code = "#ifdef _OPENMP
#               uintptr_t const oldstack = R_CStackLimit;
#               R_CStackLimit = (uintptr_t) - 1;
#             #endif
#             int const nrow    = INTEGER(GET_DIM(vDoseDropout))[0];
#             int const ncol    = INTEGER(GET_DIM(vDoseDropout))[1];
#             double const rate = *REAL(dropin);
#             int const * const condA_first  = LOGICAL(condA);
#             int const * const condB_first  = LOGICAL(condB);
#             int const * const zero_ptr     = LOGICAL(zeroAll);
#             double const * const vdose_ptr = REAL(vDoseDropout);
#             double * const out_ptr         = REAL(input);
#             double const * const freq_ptr  = REAL(freqMat);
#             
#             #pragma omp parallel for
#             for(int j=0; j < ncol; ++j) 
#             {
#               int const v = j * nrow;
#               int const *condA_ptr = condA_first;
#               int const *condB_ptr = condB_first;
#               for(int i=0; i < nrow; ++i, ++condA_ptr, ++condB_ptr) 
#               {
#                 if(*(zero_ptr + v + i)) {
#                   if(*condA_ptr)     
#                     *(out_ptr + j) *= 1.0 - (*(freq_ptr + v + i)) * rate; 
#                   else if(*condB_ptr)
#                     *(out_ptr + j) *= (*(freq_ptr + v + i)) * rate; 
#                 } else {
#                   if(*condA_ptr)      
#                     *(out_ptr + j) *= *(vdose_ptr + v + i);
#                   else if(*condB_ptr) 
#                     *(out_ptr + j) *= 1.0 - *(vdose_ptr + v + i);
#                 }
#               }
#             }
#             #ifdef _OPENMP
#                R_CStackLimit = oldstack;
#             #endif
#             
#             return R_NilValue;"
#   )
# }

probabilities.function <- function(scenario, cons, doR=FALSE) {
  # Creates a probability function.
  #
  # The probability function computes all probabilities associated with a
  # specific set of genotypes. 
  #
  # In practice, this function makes is easy to substitute C vs R
  # implementations, as well as specialize the C implementations. 

# if(require("inline") && (!scenario$doDropin) && !doR) {
#   probabilities.no.dropin.CR <- function(res, vDoseDropout, csp, unc, ...) {
#       if(length(res) != ncol(vDoseDropout))
#         stop("output vector and vDoseDropout have incompatible sizes.")
#       if(length(csp) != nrow(vDoseDropout))
#         stop("csp and vDoseDropout have incompatible sizes.")
#       if(length(unc) != nrow(vDoseDropout))
#         stop("unc and vDoseDropout have incompatible sizes.")
#       if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
#         stop("expected vDoseDropout and zeroAll to match.")
#       .probabilities.no.dropin.C(res, vDoseDropout, !csp & !unc, csp, cons$zeroAll)
#       res
#   }
#   return(probabilities.no.dropin.CR) 
# } else if(require("inline") && scenario$doDropin && !doR) {
#   probabilities.with.dropin.CR <- function(res, vDoseDropout, csp, unc, rate) {
#      if(length(res) != ncol(vDoseDropout))
#        stop("output vector and vDoseDropout have incompatible sizes.")
#      if(length(csp) != nrow(vDoseDropout))
#        stop("csp and vDoseDropout have incompatible sizes.")
#      if(length(unc) != nrow(vDoseDropout))
#        stop("unc and vDoseDropout have incompatible sizes.")
#      if(any(dim(vDoseDropout) != dim(cons$zeroAll)))
#        stop("expected vDoseDropout and zeroAll to match.")
#      if(any(dim(cons$freqMat) != dim(cons$zeroAll)))
#        stop("expected vDoseDropout and zeroAll to match.")
#      .probabilities.with.dropin.C(res, vDoseDropout, !csp & !unc, csp,
#                                   cons$zeroAll, cons$freqMat, rate)
#      res
#    }
#   return(probabilities.with.dropin.CR)
# }
  # Otherwise, return an R function. 
  function(res, vDoseDropout, csp, unc, rate) {
     res <- res                                             *
            selective.col.prod(!csp & !unc, vDoseDropout)   *
            selective.col.prod(csp, 1 - vDoseDropout)  
    if(scenario$doDropin) 
      res <- res                                                           *
             selective.col.prod(csp & cons$zeroAll, rate * cons$freqMat)   *
             selective.col.prod(!csp & !unc & cons$zeroAll,
                                1 - rate * cons$freqMat)
    res
  }
}

# if(require("inline")) {
#   .all.genotypes.per.locus.C <- compile.C(
#     signature(nContrib="integer", comb="array"),
#     code = "#ifdef _OPENMP
#               uintptr_t const oldstack = R_CStackLimit;
#               R_CStackLimit = (uintptr_t) - 1;
#             #endif
#             int const nb = *INTEGER(nContrib);
#             int const ncomb = INTEGER(GET_DIM(comb))[1];
#             int const nrow  = nb * 2;
#             int ncol  = ncomb;
#             for(int i=1; i < nb; ++i, ncol *= ncomb); 
#             SEXP result;
#             PROTECT(result = allocMatrix(INTSXP, nrow, ncol));
#             int * const ptr = INTEGER(result);
#             int * const combptr = INTEGER(comb);
#             
#             #pragma omp parallel for 
#             for(int j=0; j < ncol; ++j) 
#             {
#               int const out_index = j * nrow;
#               for(int i=nb-1, u=j; i >= 0; --i, u /= ncomb)
#               {
#                 int const in_index = (u % ncomb) << 1; 
#                 *(ptr + out_index + i * 2) = *(combptr + in_index);
#                 *(ptr + out_index + i * 2 + 1) = *(combptr + in_index + 1);
#               }
#             }
#             
#             UNPROTECT(1);
#             #ifdef _OPENMP
#                R_CStackLimit = oldstack;
#             #endif
#             return result;"
#   )
#   all.genotypes.per.locus.C <- function(nAlleles, nContrib=1) {
#     if(nAlleles < 0) stop("Negative or null number of alleles.")
#     if(nContrib < 0) stop("Negative number of contributors.")
#     if(nContrib == 0 || nAlleles == 0) return(matrix(nrow=0, ncol=0))
#     nContrib = as.integer(nContrib)
#     a = t(combinations(nAlleles, 2, repeats.allowed=TRUE))
#     if(is.null(nContrib) || is.null(a)) stop("something went very wrong.")
#     .all.genotypes.per.locus.C(nContrib, a)
#   }
#   all.genotypes.per.locus = all.genotypes.per.locus.C
# }

all.genotypes.per.locus.C <- function(nAlleles, nContrib=1) {
  if(nAlleles < 0) stop("Negative or null number of alleles.")
  if(nContrib < 0) stop("Negative number of contributors.")
  if(nContrib == 0 || nAlleles == 0) return(matrix(nrow=0, ncol=0))
  nContrib = as.integer(nContrib)
  a = t(combinations(nAlleles, 2, repeats.allowed=TRUE))
  if(is.null(nContrib) || is.null(a)) stop("something went very wrong.")
  .Call(.all.genotypes.per.locus.C, nContrib, a, PACKAGE="likeLTD")
}
all.genotypes.per.locus = all.genotypes.per.locus.C
