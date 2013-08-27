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
# ifdef _OPENMP
    double const nthreads(omp_get_max_threads()); 
# else
    double const nthreads(0);
# endif
  SEXP result;
  PROTECT(result = allocVector(REALSXP, 1));
  *static_cast<double*>(REAL(result)) = nthreads;
  UNPROTECT(1);

  return result;
}

SEXP set_nbthreads(SEXP _n)
{
# ifdef _OPENMP
    int const n = *INTEGER(_n);
    omp_set_num_threads(n); 
# endif

  return R_NilValue;
}
