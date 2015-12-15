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

    //! \brief combine doses with single stutter
	inline std::vector<genoStruct> combineDosesS(std::vector<double> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS);	                                            
	//! \brief combine doses with single and double stutter
	inline std::vector<genoStruct> combineDosesSD(std::vector<double> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd);
	//! \brief combine doses with single and over stutter
	inline std::vector<genoStruct> combineDosesSO(std::vector<double> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSo);
	//! \brief combine doses with single, double and over stutter
	inline std::vector<genoStruct> combineDosesSDO(std::vector<double> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd,std::vector<genoStruct> muSo);

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

    //! \brief modify csp
inline std::vector<cspStruct> modifyCSP(std::vector<cspStruct> csp,std::vector<double> allPosVec);

    //! \brief test cdf
    SEXP testCDF(SEXP S, SEXP Z);

    //! \brief test new cdf
    SEXP testNewCDF(SEXP S, SEXP Z);

    //! \brief test pdf
    SEXP testPDF(SEXP X, SEXP A, SEXP B);

    
    //! \brief get mean dose with x-1, x-2 and x+1 stutter
    inline std::vector<genoStruct> getDoseSDO(std::vector<double> genotypeVec, 
                        std::vector<double> stutterPosVec,std::vector<double> doubleStutterVec,
                        std::vector<double> overStutterVec,std::vector<double> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS, double meanD, 
                        double meanO,double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

    //! \brief get mean dose with x-1 and x-2 stutter
    inline std::vector<genoStruct> getDoseSD(std::vector<double> genotypeVec, 
                        std::vector<double> stutterPosVec,std::vector<double> doubleStutterVec,
                        std::vector<double> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS, double meanD, 
                        double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

       //! \brief get mean dose with x-1 and x+1 stutter
    inline std::vector<genoStruct> getDoseSO(std::vector<double> genotypeVec, 
                        std::vector<double> stutterPosVec,
                        std::vector<double> overStutterVec,std::vector<double> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS,
                        double meanO,double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

       //! \brief get mean dose with x-1 stutter
    inline std::vector<genoStruct> getDoseS(std::vector<double> genotypeVec, 
                        std::vector<double> stutterPosVec,std::vector<double> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS,
                        double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

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

//! \brief get mean dose with x-1, x-2 or x+1 stutter with possible dropin
SEXP getProbabilities(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg);


#ifdef __cplusplus
}
#endif
#endif
