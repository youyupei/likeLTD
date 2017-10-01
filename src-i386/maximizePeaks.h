#ifndef __LIKELTD_MAXIMIZEPEAKS_H_INCLUDED__  
#define __LIKELTD_MAXIMIZEPEAKS_H_INCLUDED__  
 

#include "config.h"

#include <R.h>
#include <Rdefines.h>
#include <vector>
 
#ifdef __cplusplus
extern "C" {
#endif
	
	//! \brief structure to contain genotype and associated variable
        struct genoStruct
	        {
	        double dose;
	        double genotype;
	        } ;

	//! \brief structure to contain all data associated with a CSP peak
	struct cspStruct
		{
		double heights;
		double alleles;
		double sizes;
		} ;

	//! \brief function to round a number
	inline int myRound(double d);

	//! \brief log gamma function
	inline double gammln(double xx);

	//! \brief log pdf of gamma distribution
	inline double gammalog(double x, double a, double b);

	//! \brief function to know which of gammaMu to not include
	inline bool shouldBeRemoved( genoStruct g, std::vector<cspStruct> csp);

	//! \brief log gamma function (new)
	inline double kf_lgamma(double z);

	//! \brief regularized lower incomplete gamma function, by series expansion
	double kf_gammap(double s, double z);

	//! \brief log regularized lower incomplete gamma function, by series expansion	
	double ln_kf_gammap(double s, double z);

    //! \brief get density
    inline double getDensity(std::vector<genoStruct> gammaMuVec, std::vector<cspStruct> cspModify, double scale, double cdfArg, double pdfArg);

    //! \brief test cdf
    SEXP testCDF(SEXP S, SEXP Z);

    //! \brief test new cdf
    SEXP testNewCDF(SEXP S, SEXP Z);

    //! \brief test pdf
    SEXP testPDF(SEXP X, SEXP A, SEXP B);

    
//! \brief get mean dose with x-1, x-2 and x+1 stutter with dropin
SEXP getProbabilitiesSDO_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg);
//! \brief get mean dose with x-1 and x+1 stutter with dropin
SEXP getProbabilitiesSO_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg);
//! \brief get mean dose with x-1 and x-2 stutter with dropin
SEXP getProbabilitiesSD_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg);
//! \brief get mean dose with x-1 stutter with dropin
SEXP getProbabilitiesS_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg);

//! \brief get mean dose with x-1, x-2 and x+1 stutter
SEXP getProbabilitiesSDO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals);
//! \brief get mean dose with x-1 and x+1 stutter
SEXP getProbabilitiesSO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals);
//! \brief get mean dose with x-1 and x-2 stutter
SEXP getProbabilitiesSD(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals);
//! \brief get mean dose with x-1 stutter
SEXP getProbabilitiesS(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals);


#ifdef __cplusplus
}
#endif
#endif
