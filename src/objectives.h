#ifndef LIKELTD_OBJECTIVES_H
#define LIKELTD_OBJECTIVES_H

#include <R.h>
#include <Rdefines.h>
#ifdef __cplusplus
extern "C" {
#endif
  //! Computes probabilities with dropout only.
  SEXP probabilitiesNoDropin(SEXP res, SEXP vDoseDropout, SEXP condA, 
                             SEXP condB, SEXP zeroAll);
  
  //! Computes probabilities with dropin and dropout.
  SEXP probabilitiesWithDropin(SEXP res, SEXP vDoseDropout,
                               SEXP condA, SEXP condB, SEXP zeroAll, 
                               SEXP freqMat, SEXP rate);
#ifdef __cplusplus
}
#endif
#endif
