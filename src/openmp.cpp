#include "config.h"
#include "openmp.h"

#include <R_ext/Error.h>

#ifdef _OPENMP
#  include <omp.h>
#  ifdef OPENMP_STACK
#    define CSTACK_DEFNS 7
#    include <Rinterface.h>
#  endif
#endif

SEXP nbthreads()
{
# ifdef OPENMP_STACK
//    uintptr_t const oldstack = R_CStackLimit;
//    R_CStackLimit = (uintptr_t) - 1;
# endif
  SEXP result;
  PROTECT(result = allocVector(REALSXP, 1));
  double *__result = REAL(result);
# ifdef _OPENMP
    *__result = omp_get_max_threads(); 
# else
    *__result = 0;
# endif

# ifdef OPENMP_STACK
//    R_CStackLimit = oldstack;
# endif
  UNPROTECT(1);
  return result;
}
