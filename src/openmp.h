#ifndef LIKELTD_OPENMP_H
#define LIKELTD_OPENMP_H

#include "config.h"

#include <R.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif
  //! \brief Returns number of treads
  SEXP nbthreads();
#ifdef __cplusplus
}
#endif
#endif
