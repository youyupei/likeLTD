#include "config.h"
#include <R_ext/Rdynload.h>
#include "genetics.h"
#include "objectives.h"
#include "openmp.h"
#include "maximizePeaks.h"

R_CallMethodDef callMethods[]  = {
       {"allGenotypesPerLocus", (DL_FUNC) &allGenotypesPerLocus, 2},
       {"addProfilesToEPG", (DL_FUNC) &addProfilesToEPG, 3},
       {"probabilitiesWithDropin", (DL_FUNC) &probabilitiesWithDropin, 7},
       {"probabilitiesNoDropin", (DL_FUNC) &probabilitiesNoDropin, 5},
       {"powerAdjustment", (DL_FUNC) &powerAdjustment, 4},
       {"doseFraction", (DL_FUNC) &doseFraction, 3},
       {"emptyAlleles", (DL_FUNC) &emptyAlleles, 2},
       {"fractionsAndHet", (DL_FUNC) &fractionsAndHet, 2},
       {"relatednessFactors", (DL_FUNC) &relatednessFactors, 5},
       {"nbthreads", (DL_FUNC) &nbthreads, 0},
       {"set_nbthreads", (DL_FUNC) &set_nbthreads, 1},
		{"getProbabilitiesSDO", (DL_FUNC) &getProbabilitiesSDO, 16},
		{"getProbabilitiesSD", (DL_FUNC) &getProbabilitiesSDO, 15},
		{"getProbabilitiesSO", (DL_FUNC) &getProbabilitiesSDO, 15},
		{"getProbabilitiesS", (DL_FUNC) &getProbabilitiesSDO, 14},
		{"getProbabilitiesSDO_dropin", (DL_FUNC) &getProbabilitiesSDO_dropin, 18},
		//{"getProbabilitiesSDO_dropinOptimised", (DL_FUNC) &getProbabilitiesSDO_dropinOptimised, 18},
//		{"testCDF", (DL_FUNC) &testCDF, 2},
//		{"testPDF", (DL_FUNC) &testPDF, 3},
       {NULL, NULL, 0},
       {NULL, NULL, 0},
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

