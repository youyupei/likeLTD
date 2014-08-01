//==================================
// define class functions
 
//================================
// include guard
#ifndef __LIKELTD_MAXIMIZEPEAKS_h_INCLUDED__  
#define __LIKELTD_MAXIMIZEPEAKS_H_INCLUDED__  
 
//================================
// forward declared dependencies
 
//=================================
// included dependencies

#include "config.h"

#include <R.h>
#include <Rdefines.h>


 
#ifdef __cplusplus
extern "C" {
#endif
	//! \brief peak mean dose
        SEXP peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames);
  
#ifdef __cplusplus
}
#endif
#endif
