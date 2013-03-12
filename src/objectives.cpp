#include "objectives.h"

#include <cmath>
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

inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

// should be much more precise with large b
inline double fastPrecisePow(double a, double b) {
  // calculate approximation with fraction of the exponent
  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
 
  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }
 
  return r * u.d;
}

// Computes adjustment + exponential. 
SEXP tvedebrinkAdjustment(SEXP allEPG, SEXP zeroAll, SEXP localAdjustment, 
                          SEXP tvedebrink)
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(allEPG))[0];
  int const ncol = INTEGER(GET_DIM(allEPG))[1];
  int const matsize = nrow * ncol;
  if(nrow != INTEGER(GET_DIM(zeroAll))[0]) {
    error("First dimension of allEPG and zeroAll do not match.");
    return R_NilValue;
  }
  if(ncol != INTEGER(GET_DIM(zeroAll))[1]) {
    error("Second dimension of allEPG and zeroAll do not match.");
    return R_NilValue;
  }
  double      const adjust   = *REAL(localAdjustment);
  double      const tvede    = *REAL(tvedebrink);
  int const * const zero_ptr = LOGICAL(zeroAll);
  double    * const epg_ptr  = REAL(allEPG);

# pragma omp parallel for 
  for(int j=0; j < matsize; ++j) 
    if(not *(zero_ptr +j))
      *(epg_ptr + j) =  std::exp(tvede * std::log(adjust * (*(epg_ptr + j))));
// using pow is slower somehow...
//     *(epg_ptr + j) =  std::pow(adjust * (*(epg_ptr + j)), tvede);

# ifdef _OPENMP
    R_CStackLimit = oldstack;
# endif

  return R_NilValue;
}


// Computes faction allEPG * dropout / (allEPG + 1 - dropout)
SEXP fraction(SEXP allEPG, SEXP zeroAll, SEXP dropout)
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(allEPG))[0];
  int const ncol = INTEGER(GET_DIM(allEPG))[1];
  if(nrow != INTEGER(GET_DIM(zeroAll))[0]) {
    error("First dimension of allEPG and zeroAll do not match.");
    return R_NilValue;
  }
  if(ncol != INTEGER(GET_DIM(zeroAll))[1]) {
    error("Second dimension of allEPG and zeroAll do not match.");
    return R_NilValue;
  }

  SEXP result;
  PROTECT(result = allocMatrix(REALSXP, nrow, ncol));
  double         const rate     = *REAL(dropout);
  int const *    const zero_ptr = LOGICAL(zeroAll);
  double const * const epg_ptr  = REAL(allEPG);
  double       * const out_ptr  = REAL(result);

# pragma omp parallel for 
  for(int j=0; j < ncol*nrow; ++j) 
    if(not *(zero_ptr +j))
    {
      double const newfact = *(epg_ptr + j) * rate;
      *(out_ptr + j) = newfact / (newfact + 1.0 - rate);
    }

# ifdef _OPENMP
    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);

  return result;
}

