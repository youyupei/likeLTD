#ifndef LIKELTD_OBJECTIVES_H
#define LIKELTD_OBJECTIVES_H

#include "config.h"

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
  //! \details computes (locusAdjustment * allEPG[!zeroAll])^power. allEPG
  //!          is both input and output.
  SEXP powerAdjustment(SEXP allEPG, SEXP zeroAll, SEXP locusAdjustment, 
                            SEXP power);

  //! \brief Computes faction allEPG * dropout / (allEPG + 1 - dropout)
  //! \details Does the above operation only for zeroAll == false. 
  //!          Returns the result of the operation in new matrix.
  SEXP doseFraction(SEXP allEPG, SEXP zeroAll, SEXP dropout);

  //! \brief computes matrix indicating allele presence. 
  //! \details returns allZeros
  SEXP emptyAlleles(SEXP genotypes, SEXP knownZero);
 
  //! \brief Computes allele and heterozygote factors.
  //! \details The value of genotypes should never exceed the size of
  //! fractions, otherwise we'll get a segfault. To expensive to check though. 
  SEXP fractionsAndHet(SEXP genotypes, SEXP fractions);

  //! \brief factors in relatedness.
  //! \details There are only four possible factors corresponding to four
  //!          possible hypothesis: 
  //!            
  //!            (i) neither queried alleles are present.
  //!            (ii) first queried allele is in first unknown contributor, 
  //!            (iii) second queried allele is in first unknown contributor,
  //!            (iv) both queried alleles are in first unknown contributor,
  //!
  //!          These four factors are determined for and applied to each entry
  //!          in inout.
  //! \return NULL value. inout is both input and output.
  SEXP relatednessFactors(SEXP inout, SEXP relatednessR, SEXP genotypes,
                          SEXP queriedR, SEXP frequenciesR);
#ifdef __cplusplus
}
#endif
#endif
