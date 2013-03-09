#include "genetics.h"

#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  define CSTACK_DEFNS 7
#  include <Rinterface.h>
#endif

SEXP allGenotypesPerLocus(SEXP nContrib, SEXP comb)
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  int const nb = *INTEGER(nContrib);
  int const ncomb = INTEGER(GET_DIM(comb))[1];
  int const nrow  = nb * 2;
  int ncol  = ncomb;
  for(int i=1; i < nb; ++i, ncol *= ncomb);
  SEXP result;
  PROTECT(result = allocMatrix(INTSXP, nrow, ncol));
  int * const ptr = INTEGER(result);
  int * const combptr = INTEGER(comb);

# pragma omp parallel private(ncol) for
  for(int j=0; j < ncol; ++j)
  {
    int const out_index = j * nrow;
    for(int i=nb-1, u=j; i >= 0; --i, u /= ncomb)
    {
      int const in_index = (u % ncomb) << 1;
      *(ptr + out_index + i * 2) = *(combptr + in_index);
      *(ptr + out_index + i * 2 + 1) = *(combptr + in_index + 1);
    }
  }

  UNPROTECT(1);
# ifdef _OPENMP
     R_CStackLimit = oldstack;
# endif
  return result;
}
