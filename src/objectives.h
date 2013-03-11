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
  //! \brief Computes adjustment + exponential. 
  //! \details computes (localAdjustment * allEPG[!zeroAll])^tvedebrink. allEPG
  //!          is both input and output.
  SEXP tvedebrinkAdjustment(SEXP allEPG, SEXP zeroAll, SEXP localAdjustment, 
                            SEXP tvedebrink);

  // \brief Computes faction allEPG * dropout / (allEPG + 1 - dropout)
  // \details Does the above operation only for zeroAll == false. 
  //          Returns the result of the operation in new matrix.
  SEXP fraction(SEXP allEPG, SEXP zeroAll, SEXP dropout);
#ifdef __cplusplus
}
#endif
#endif
