#ifndef __LIKELTD_MAXIMIZEPEAKS_h_INCLUDED__  
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
	
	//! \brief peak mean dose
        //SEXP peakMeanDose(SEXP genotype, SEXP alleles,SEXP heights,SEXP sizes,SEXP DNAcont,SEXP stutter,SEXP degradation,SEXP fragL, SEXP fragN,SEXP repAdjust);
	std::vector<genoStruct> peakMeanDose(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec,double stutterTrue, double stutterFalse,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, int nGen, int nCSP, int nCont, int nFrag);

	//! \brief combine doses
	std::vector<genoStruct> combineDoses(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS);

    	//! \brief probability of peaks
        SEXP probabilityPeaks(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP DNAcont, SEXP stutterMean, SEXP stutterAdjust, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust, SEXP detectionThresh);

	//! \brief probability of a single genotype
	double singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutterTrue, double stutterFalse, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, int nFrag, double cdfArg, double pdfArg);
	//std::vector<double> singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, double detectionThresh, int currentComb, int nGen, int nCSP, int nCont, int nFrag);

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
double getDensity(std::vector<genoStruct> gammaMuVec, std::vector<cspStruct> cspModify, double scale, double cdfArg, double pdfArg);

//! \brief modify csp
std::vector<cspStruct> modifyCSP(std::vector<cspStruct> csp,std::vector<float> allPosVec);

#ifdef __cplusplus
}
#endif
#endif
