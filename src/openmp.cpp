#include "openmp.h"

#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  define CSTACK_DEFNS 7
#  include <Rinterface.h>
#endif

SEXP nbthreads()
{
# ifdef _OPENMP
    uintptr_t const oldstack = R_CStackLimit;
    R_CStackLimit = (uintptr_t) - 1;
# endif
  SEXP result;
  PROTECT(result = allocVector(REALSXP, 1));
  double *__result = REAL(result);
# ifdef _OPENMP
    *__result = omp_get_max_threads(); 
# else
    *__result = 0;
# endif

# ifdef _OPENMP
    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);
  return result;
}
