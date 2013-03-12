#include <R_ext/Rdynload.h>
#include "genetics.h"
#include "objectives.h"
#include "openmp.h"

R_CallMethodDef callMethods[]  = {
       {"allGenotypesPerLocus", (DL_FUNC) &allGenotypesPerLocus, 2},
       {"addProfilesToEPG", (DL_FUNC) &addProfilesToEPG, 3},
       {"probabilitiesWithDropin", (DL_FUNC) &probabilitiesWithDropin, 7},
       {"probabilitiesNoDropin", (DL_FUNC) &probabilitiesNoDropin, 5},
       {"tvedebrinkAdjustment", (DL_FUNC) &tvedebrinkAdjustment, 4},
       {"fraction", (DL_FUNC) &fraction, 3},
       {"emptyAlleles", (DL_FUNC) &emptyAlleles, 2},
       {"nbthreads", (DL_FUNC) &nbthreads, 0},
       {NULL, NULL, 0}
};

#ifdef __cplusplus
extern "C" 
{ 
#endif
  void R_init_likeLTD(DllInfo *info)
  {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  }
#ifdef __cplusplus
}
#endif
