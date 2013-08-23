#include "config.h"
#include "objectives.h"

#include <cmath>
#include <string>
#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  ifdef OPENMP_STACK
#    define CSTACK_DEFNS 7
#    include <Rinterface.h>
#  endif
#endif

// Computes probabilities with dropout only.
SEXP probabilitiesNoDropin(SEXP input, SEXP vDoseDropout, SEXP condA, SEXP condB,
                           SEXP zeroAll)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
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
      if(zero_ptr[v + i]) continue;
      if(*condA_ptr)      out_ptr[j] *= vdose_ptr[v + i];
      else if(*condB_ptr) out_ptr[j] *= 1.0 - vdose_ptr[v + i];
    }
  }
# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif

  return input;
}

// Computes probabilities with dropin and dropout.
SEXP probabilitiesWithDropin(SEXP input, SEXP vDoseDropout, SEXP condA,
                             SEXP condB, SEXP zeroAll, SEXP freqMat,
                             SEXP dropin)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
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
      if(zero_ptr[v + i]) {
        if(*condA_ptr)      out_ptr[j] *= 1.0 - freq_ptr[i] * rate; 
        else if(*condB_ptr) out_ptr[j] *= freq_ptr[i] * rate; 
      } else {
        if(*condA_ptr)      out_ptr[j] *= vdose_ptr[v + i];
        else if(*condB_ptr) out_ptr[j] *= 1.0 - vdose_ptr[v + i];
      }
    }
  }
# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif

  return input;
}

// Computes adjustment + exponential. 
SEXP powerAdjustment(SEXP allEPG, SEXP zeroAll, SEXP locusAdjustment, 
                          SEXP power)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
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
  double      const adjust   = *REAL(locusAdjustment);
  double      const tvede    = *REAL(power);
  int const * const zero_ptr = LOGICAL(zeroAll);
  double    * const epg_ptr  = REAL(allEPG);

# pragma omp parallel for 
  for(int j=0; j < matsize; ++j) 
    if(not zero_ptr[j])
      epg_ptr[j] =  std::exp(tvede * std::log(adjust * (epg_ptr[j])));
// using pow is slower somehow...
//     epg_ptr[j] =  std::pow(adjust * (epg_ptr[j]), tvede);

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif

  return R_NilValue;
}


// Computes faction allEPG * dropout / (allEPG + 1 - dropout)
SEXP doseFraction(SEXP allEPG, SEXP zeroAll, SEXP dropout)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
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
    if(not zero_ptr[j])
    {
      double const newfact = epg_ptr[j] * rate;
      out_ptr[j] = newfact / (newfact + 1.0 - rate);
    }

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);

  return result;
}

// Computes matrix indicating presence/absence of alleles.
SEXP emptyAlleles(SEXP genotypes, SEXP knownZero)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(genotypes))[0];
  int const ncol = INTEGER(GET_DIM(genotypes))[1];
  int const nAlleles = LENGTH(knownZero);
  if(nrow % 2 != 0) {
    error("Expected first dimension of genotypes to be even.");
    return R_NilValue;
  }

  SEXP result;
  PROTECT(result = allocMatrix(LGLSXP, nAlleles, ncol));
  int const * const geno_ptr  = INTEGER(genotypes);
  int const * const known_ptr = LOGICAL(knownZero);
  int       * const out_ptr   = LOGICAL(result);
  int const cpysize = sizeof(int) * nAlleles;

# pragma omp for
  for(int j=0; j < ncol; ++j)
  {
    int const v = j * nAlleles;
    // copies zeros from known profiles to thos column
    memcpy((void*) (out_ptr + v), (void*) known_ptr, cpysize);
    // now sets additional alleles to zero
    for(int i=j*nrow; i < (j + 1) * nrow; ++i)
      out_ptr[v + geno_ptr[i] - 1] = 0;
  }

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);

  return result;
}

//! Computes allele fractions.
SEXP fractionsAndHet(SEXP genotypes, SEXP fractions)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(genotypes))[0];
  int const ncol = INTEGER(GET_DIM(genotypes))[1];
  if(nrow % 2 != 0) {
    error("Expected first dimension of genotypes to be even.");
    return R_NilValue;
  }

  SEXP result;
  PROTECT(result = allocVector(REALSXP, ncol));

  int const * const geno_ptr    = INTEGER(genotypes);
  double const * const frac_ptr = REAL(fractions);
  double       * const out_ptr  = REAL(result);

# pragma omp for
  for(int j=0; j < ncol; ++j)
  {
    int const * i_geno = geno_ptr + j*nrow;
    // compute heterozygote factor.
    int n = 1;
    for(int i=0; i < nrow; i += 2, i_geno += 2)
      if(*i_geno < *(i_geno + 1)) n *= 2; 

    // compute allele fraction factor.
    i_geno = geno_ptr + j*nrow;
    out_ptr[j] = double(n);
    for(int i=0; i < nrow; ++i, ++i_geno)
      out_ptr[j] *= frac_ptr[*i_geno - 1];
  }

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);

  return result;
}

// Computes relatedness factors and multiplies to inout.
SEXP relatednessFactors(SEXP inout, SEXP relatednessR, SEXP genotypes, 
                        SEXP queriedR, SEXP frequenciesR)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow = INTEGER(GET_DIM(genotypes))[0];
  int const ncol = INTEGER(GET_DIM(genotypes))[1];
  if(nrow % 2 != 0) {
    error("Expected first dimension of genotypes to be even.");
    return R_NilValue;
  }
  if(nrow < 2) {
    error("Expected first dimension of genotypes to be at least 2.");
    return R_NilValue;
  }
  if(ncol != LENGTH(inout)) {
    error("Expected first dimension of genotypes and length of inout to "
          "match.");
    return R_NilValue;
  }
  if(LENGTH(queriedR) != 2) {
    error("Expected queried to be an integer with two elements.");
    return R_NilValue;
  }
  if(LENGTH(relatednessR) != 2) {
    error("Expected relatednessR to be a double with two elements.");
    return R_NilValue;
  }
  if(LENGTH(frequenciesR) != 2) {
    error("Expected relatednessR to be a double with two elements.");
    return R_NilValue;
  }

  // Translate input constants to C.
  double const relatedness[2] = { REAL(relatednessR)[0],   
                                  REAL(relatednessR)[1]  }; 
  double const frequencies[2] = { REAL(frequenciesR)[0],   
                                  REAL(frequenciesR)[1]  }; 
  int const queried[2] = { INTEGER(queriedR)[0], INTEGER(queriedR)[1]  }; 

  // The following builds up the factors which will modify the inout vector. 
  // There are only four possible factors corresponding to four possible
  // hypothesis: (i) neither queried alleles are present.
  //            (ii) first queried allele is in first unknown contributor, 
  //            (iii) second queried allele is in first unknown contributor,
  //            (iv) both queried alleles are in first unknown contributor,
  // The final part of  the factors depends on whether that particular genotype
  // is homozygote or not. This is checked within the loop.
  double const factor0 = 1e0 - relatedness[0] - relatedness[1];
  double const factors[3] = { 0.5 * relatedness[0] / frequencies[0],
                              0.5 * relatedness[0] / frequencies[1], 
                              relatedness[1] / frequencies[0]
                                             / frequencies[1] };
  // define pointers to input/output arrays.
  int const * const geno_ptr    = INTEGER(genotypes);
  double       * const out_ptr  = REAL(inout);

  // Finally loop.
# pragma omp for
  for(int j=0; j < ncol; ++j)
  {
    int const *i_geno = geno_ptr + j*nrow; 

    bool const first  = i_geno[0] == queried[0] or i_geno[1] == queried[0];
    bool const second = i_geno[0] == queried[1] or i_geno[1] == queried[1];
    bool const is_homozygote = i_geno[0] == i_geno[1];

    double result = factor0;
    if(first)  result += is_homozygote ? factors[0]: factors[0] * 0.5;
    if(second) result += is_homozygote ? factors[1]: factors[1] * 0.5;
    if(first and second)
      result += is_homozygote ? factors[2]: factors[2] * 0.5;

    out_ptr[j] *= result;
  }

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif

  return R_NilValue;
}
