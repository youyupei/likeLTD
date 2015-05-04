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
	        float dose;
	        float genotype;
	        } ;

	//! \brief structure to contain all data associated with a CSP peak
	struct cspStruct
		{
		float heights;
		float alleles;
		float sizes;
		} ;

    //! \brief combine doses with single stutter
	inline std::vector<genoStruct> combineDosesS(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS);	                                            
	//! \brief combine doses with single and double stutter
	inline std::vector<genoStruct> combineDosesSD(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd);
	//! \brief combine doses with single and over stutter
	inline std::vector<genoStruct> combineDosesSO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSo);
	//! \brief combine doses with single, double and over stutter
	inline std::vector<genoStruct> combineDosesSDO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd,std::vector<genoStruct> muSo);

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
	static double kf_gammap(double s, double z);

	//! \brief log regularized lower incomplete gamma function, by series expansion	
	static double ln_kf_gammap(double s, double z);

    //! \brief get density
    inline double getDensity(std::vector<genoStruct> gammaMuVec, std::vector<cspStruct> cspModify, double scale, double cdfArg, double pdfArg);

    //! \brief modify csp
inline std::vector<cspStruct> modifyCSP(std::vector<cspStruct> csp,std::vector<float> allPosVec);

    //! \brief test cdf
    SEXP testCDF(SEXP S, SEXP Z);

    //! \brief test pdf
    SEXP testPDF(SEXP X, SEXP A, SEXP B);

    
    //! \brief get mean dose with x-1, x-2 and x+1 stutter
    inline std::vector<genoStruct> getDoseSDO(std::vector<float> genotypeVec, 
                        std::vector<float> stutterPosVec,std::vector<float> doubleStutterVec,
                        std::vector<float> overStutterVec,std::vector<float> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS, double meanD, 
                        double meanO,double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

    //! \brief get mean dose with x-1 and x-2 stutter
    inline std::vector<genoStruct> getDoseSD(std::vector<float> genotypeVec, 
                        std::vector<float> stutterPosVec,std::vector<float> doubleStutterVec,
                        std::vector<float> overStutterVec,std::vector<float> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS, double meanD, 
                        double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

       //! \brief get mean dose with x-1 and x+1 stutter
    inline std::vector<genoStruct> getDoseSO(std::vector<float> genotypeVec, 
                        std::vector<float> stutterPosVec,std::vector<float> doubleStutterVec,
                        std::vector<float> overStutterVec,std::vector<float> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS,
                        double meanO,double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

       //! \brief get mean dose with x-1 stutter
    inline std::vector<genoStruct> getDoseS(std::vector<float> genotypeVec, 
                        std::vector<float> stutterPosVec,std::vector<float> doubleStutterVec,
                        std::vector<float> overStutterVec,std::vector<float> allPosVec, 
                        std::vector<double> DNAcontVec, double gradientS,
                        double interceptS, std::vector<double> degVec,
                        std::vector<double> fragVecL, std::vector<double> fragVecN, 
                        std::vector<double> stutterIndex, int nGen, int nFrag);

    //! get probabilities of genotype combinations with x-1, x-2 and x+1 stutter
    SEXP getProbabilitiesSDO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, 
                    SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, 
                    SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, 
                    SEXP detectionThresh, SEXP databaseVals);

    //! get probabilities of genotype combinations with x-1 and x-2 stutter
    SEXP getProbabilitiesSD(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, 
                    SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, 
                    SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, 
                    SEXP detectionThresh, SEXP databaseVals);

    //! get probabilities of genotype combinations with x-1 and x+1 stutter
    SEXP getProbabilitiesSO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanO, 
                    SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, 
                    SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, 
                    SEXP detectionThresh, SEXP databaseVals);

    //! get probabilities of genotype combinations with x-1 stutter
    SEXP getProbabilitiesS(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP interceptS, 
                    SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, 
                    SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, 
                    SEXP detectionThresh, SEXP databaseVals);

    //! get probabilities of genotype combinations with x-1, x-2 and x+1 stutter allowing for dropin
SEXP getProbabilitiesSDO_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP meanD, SEXP meanO, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP stutterVals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin);

    //! get probabilities of genotype combinations with x-1, x-2 and x+1 stutter allowing for dropin
//SEXP getProbabilitiesSDO_dropinOptimised(SEXP genotypeArray, SEXP DNAcont, SEXP meanD, SEXP meanO, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP stutterVals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin);


#ifdef __cplusplus
}
#endif
#endif
