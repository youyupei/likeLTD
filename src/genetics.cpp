#include "config.h"
#include "genetics.h"

#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  ifdef OPENMP_STACK
#    define CSTACK_DEFNS 7
//#    include <Rinterface.h>
#  endif
#endif

SEXP allGenotypesPerLocus(SEXP nContrib, SEXP comb)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nb = *INTEGER(nContrib);
  int const ncomb = INTEGER(GET_DIM(comb))[1];
  int const nrow  = nb * 2;
  int nncol  = ncomb;
  for(int i=1; i < nb; ++i, nncol *= ncomb);
  int const ncol = nncol;
  SEXP result;
  PROTECT(result = allocMatrix(INTSXP, nrow, ncol));
  if(nrow and ncol)
  {
    int * const ptr = INTEGER(result);
    int * const combptr = INTEGER(comb);
    
#   pragma omp parallel for
    for(int j=0; j < ncol; ++j)
    {
      int const out_index = j * nrow;
      for(int i=nb-1, u=j; i >= 0; --i, u /= ncomb)
      {
        int const l = (u % ncomb) << 1;
        int const m = out_index + (i << 1);
        ptr[m] = combptr[l];
        ptr[m + 1] = combptr[l + 1];
      }
    }
  }
# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);
  return result;
}

SEXP addProfilesToEPG(SEXP allEPG, SEXP genotypes, SEXP doses)
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nrow      = INTEGER(GET_DIM(allEPG))[0];
  int const ncol      = INTEGER(GET_DIM(allEPG))[1];
  int const nUnknowns = INTEGER(GET_DIM(genotypes))[0];
  if(ncol != INTEGER(GET_DIM(genotypes))[1]) 
  {
    error("Second dimension of allEPG and genotypes do not match.");
    return R_NilValue;
  }
  if(INTEGER(GET_DIM(doses))[1] != nUnknowns) 
  {
    error("Second dimension of doses does not match genotypes.");
    return R_NilValue;
  }
  if(nrow != INTEGER(GET_DIM(doses))[0]) 
  {
    error("First dimension of allEPG and doses do not match.");
    return R_NilValue;
  }

  double       * const epg_ptr    = REAL(allEPG);
  int const    * const geno_first = INTEGER(genotypes);
  double const * const doses_ptr  = REAL(doses);
  if(nrow and ncol)
  {
#   pragma omp parallel for
    for(int j=0; j < ncol; ++j)
    {
      int const column = j * nrow;
      int const *geno_ptr = geno_first + j * nUnknowns;
      for(int u=0, v=0; u < nUnknowns; ++u, v += nrow, ++geno_ptr) 
      {
        int const index = *geno_ptr - 1;
        epg_ptr[column + index] += doses_ptr[index + v];
      } // loop over unknown contributors.
    } // loop over space of compatible genotypes
  }
# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  return R_NilValue;
}
