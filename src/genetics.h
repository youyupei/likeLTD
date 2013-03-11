#ifndef LIKELTD_GENETICS_H
#define LIKELTD_GENETICS_H

#include <R.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif
  //! \brief Computes all possible genotypes per locus.
  //! \param[in] nContrib: Number of contributors.
  //! \param[in] comb: All possible genotypes for a single unknown
  //!                  contributor.
  SEXP allGenotypesPerLocus(SEXP nContrib, SEXP comb);


  //! \brief Helper function to speedup compute of allEPG.
  //! \param[inout] allEPG: input/output matrix to which doses from unknown
  //!               contributors are added.
  //! \param[in] genotypes: All possible genotype combinations from unknown
  //!                       contributors.
  //! \param[in] doses: Matrix with doses for each allele (rows) and each
  //!                   contributor (columns), e.g. 
  //!                   doses_i,j = r_j * (1+d_j)^f_i.
  //! \details Called in all.epg.per.locus.
  //!          In R code, it would do the following at place of call:
  //!
  //!          for(j in 1:ncol(genotypes)) {
  //!            for(u in 1:nUnknowns) {
  //!              index = genotypes[u, j]
  //!              allEPG[index, j] = allEPG[index, j]  + unknownDoses[index, u]
  //!            }
  //!          }
  SEXP addProfilesToEPG(SEXP allEPG, SEXP genotypes, SEXP fragLengths);
#ifdef __cplusplus
}
#endif
#endif
