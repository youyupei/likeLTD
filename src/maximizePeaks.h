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
#include <vector>




 
#ifdef __cplusplus
extern "C" {
#endif

        struct genoStruct
	        {
	        double genotype;
	        double dose;
	        } ;

	struct cspStruct
		{
		double alleles;
		double heights;
		double sizes;
		} ;
	
	//! \brief peak mean dose
        SEXP peakMeanDose(SEXP genotype, SEXP alleles,SEXP heights,SEXP sizes,SEXP DNAcont,SEXP stutter,SEXP degradation,SEXP fragL, SEXP fragN,SEXP repAdjust);

    //! \brief probability of peaks
    //    SEXP probabilityPeaksCPP(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP sizes, SEXP DNAcont, SEXP stutter, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust);

    //! \brief density of gamma distribution
        double gammaDensity(double x,double k,double theta);
	
	//! \brief gamma function
	double gamm(double x);

	// brief pdf of gamma distribution
	double gammaPdf(double x, double a, double b);

    //! \brief function to know which of gammaMu to not include
	bool shouldBeRemoved( genoStruct g, std::vector<cspStruct> csp);
#ifdef __cplusplus
}
#endif
#endif
