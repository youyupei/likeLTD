#include "objectives.h"

#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  define CSTACK_DEFNS 7
#  include <Rinterface.h>
#endif

// Computes probabilities with dropout only.
SEXP probabilitiesNoDropin(SEXP input, SEXP vDoseDropout, SEXP condA, SEXP condB,
                           SEXP zeroAll)
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(vDoseDropout))[0];
  int const ncol = INTEGER(GET_DIM(vDoseDropout))[1];
  int const * const condA_first  = LOGICAL(condA);
  int const * const condB_first  = LOGICAL(condB);
  int const * const zero_ptr     = LOGICAL(zeroAll);
  double const * const vdose_ptr = REAL(vDoseDropout);
  double * const out_ptr         = REAL(input);

# pragma omp parallel for 
  for(int j=0; j < ncol; ++j) 
  {
    int const v = j * nrow;
    int const *condA_ptr = condA_first;
    int const *condB_ptr = condB_first;
    for(int i=0; i < nrow; ++i, ++condA_ptr, ++condB_ptr) 
    {
      if(*(zero_ptr + v + i)) continue;
      if(*condA_ptr)      
        *(out_ptr + j) *= *(vdose_ptr + v + i);
      else if(*condB_ptr) 
        *(out_ptr + j) *= 1.0 - *(vdose_ptr + v + i);
    }
  }
# ifdef _OPENMP
    R_CStackLimit = oldstack;
# endif

  return input;
}

// Computes probabilities with dropin and dropout.
SEXP probabilitiesWithDropin(SEXP input, SEXP vDoseDropout, SEXP condA,
                             SEXP condB, SEXP zeroAll, SEXP freqMat,
                             SEXP dropin)
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow    = INTEGER(GET_DIM(vDoseDropout))[0];
  int const ncol    = INTEGER(GET_DIM(vDoseDropout))[1];
  double const rate = *REAL(dropin);
  int const * const condA_first  = LOGICAL(condA);
  int const * const condB_first  = LOGICAL(condB);
  int const * const zero_ptr     = LOGICAL(zeroAll);
  double const * const vdose_ptr = REAL(vDoseDropout);
  double * const out_ptr         = REAL(input);
  double const * const freq_ptr  = REAL(freqMat);
  
# pragma omp parallel for
  for(int j=0; j < ncol; ++j) 
  {
    int const v = j * nrow;
    int const *condA_ptr = condA_first;
    int const *condB_ptr = condB_first;
    for(int i=0; i < nrow; ++i, ++condA_ptr, ++condB_ptr) 
    {
      if(*(zero_ptr + v + i)) {
        if(*condA_ptr)     
          *(out_ptr + j) *= 1.0 - (*(freq_ptr + i)) * rate; 
        else if(*condB_ptr)
          *(out_ptr + j) *= (*(freq_ptr + i)) * rate; 
      } else {
        if(*condA_ptr)      
          *(out_ptr + j) *= *(vdose_ptr + v + i);
        else if(*condB_ptr) 
          *(out_ptr + j) *= 1.0 - *(vdose_ptr + v + i);
      }
    }
  }
# ifdef _OPENMP
    R_CStackLimit = oldstack;
# endif

  return input;
}
