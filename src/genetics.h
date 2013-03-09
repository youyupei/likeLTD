#ifndef LIKELTD_GENETICS_H
#define LIKELTD_GENETICS_H

#include <R.h>
#include <Rdefines.h>
extern "C" {
  //! \brief Computes all possible genotypes per locus.
  //! \param[in] nContrib: Number of contributors.
  //! \param[in] comb: All possible genotypes for a single unknown
  //!                  contributor.
  SEXP allGenotypesPerLocus(SEXP nContrib, SEXP comb);
}
#endif
