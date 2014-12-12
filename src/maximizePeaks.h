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
	//! \brief peak mean dose without single stutter
	inline std::vector<genoStruct> peakMeanDoseS(std::vector<float> genotypeVec, std::vector<float> stutterPosVec,  
	                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, 
	                                            std::vector<double> DNAcontVec,double stutterMean, double stutterAdjust, 
	                                            double stutterGradient,
	                                            std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
	                                            int nGen, int nCSP, int nCont, int nFrag);
	//! \brief peak mean dose with single and double stutter
	inline std::vector<genoStruct> peakMeanDoseSD(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, 
	                                            std::vector<float> doubleStutterVec, 
	                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, 
	                                            std::vector<double> DNAcontVec,double stutterMean, double stutterAdjust, 
	                                            double doubleStutterRate,
	                                            double stutterGradient,
	                                            std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
	                                            int nGen, int nCSP, int nCont, int nFrag);
	//! \brief peak mean dose with single and over stutter
	inline std::vector<genoStruct> peakMeanDoseSO(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, 
	                                            std::vector<float> overStutterVec, 
	                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, 
	                                            std::vector<double> DNAcontVec,double stutterMean, double stutterAdjust, 
	                                            double overStutterRate,
	                                            double stutterGradient,
	                                            std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
	                                            int nGen, int nCSP, int nCont, int nFrag);
	//! \brief peak mean dose with single, double and over stutter
	inline std::vector<genoStruct> peakMeanDoseSDO(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, 
	                                            std::vector<float> doubleStutterVec, std::vector<float> overStutterVec, 
	                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, 
	                                            std::vector<double> DNAcontVec,double stutterMean, double stutterAdjust, 
	                                            double doubleStutterRate, double overStutterRate,
	                                            double stutterGradient,
	                                            std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
	                                            int nGen, int nCSP, int nCont, int nFrag);

    //! \brief combine doses with single stutter
	inline std::vector<genoStruct> combineDosesS(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS);	                                            
	//! \brief combine doses with single and double stutter
	inline std::vector<genoStruct> combineDosesSD(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd);
	//! \brief combine doses with single and over stutter
	inline std::vector<genoStruct> combineDosesSO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSo);
	//! \brief combine doses with single, double and over stutter
	inline std::vector<genoStruct> combineDosesSDO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS,std::vector<genoStruct> muSd,std::vector<genoStruct> muSo);

	//! \brief probability of a single genotype without single stutter
	inline double singleGenotypeS(std::vector<double> genotypeArray, std::vector<cspStruct> csp, 
	                            std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double stutterGradient, 
	                            double repAdjust,double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
	                            std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, 
	                            int nFrag, double cdfArg, double pdfArg);	  
	//! \brief probability of a single genotype with single and double stutter
	inline double singleGenotypeSD(std::vector<double> genotypeArray, std::vector<cspStruct> csp, 
	                            std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double doubleStutterRate, double stutterGradient, 
	                            double repAdjust,double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
	                            std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, 
	                            int nFrag, double cdfArg, double pdfArg);
	//! \brief probability of a single genotype with single and over stutter
	inline double singleGenotypeSO(std::vector<double> genotypeArray, std::vector<cspStruct> csp, 
	                            std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double overStutterRate, double stutterGradient,
	                            double repAdjust,double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
	                            std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, 
	                            int nFrag, double cdfArg, double pdfArg);
	//! \brief probability of a single genotype with single, double and over stutter
	inline double singleGenotypeSDO(std::vector<double> genotypeArray, std::vector<cspStruct> csp, 
	                            std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double doubleStutterRate, double overStutterRate, double stutterGradient, 
	                            double repAdjust,double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
	                            std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, 
	                            int nFrag, double cdfArg, double pdfArg);

    
    //! \brief probability of peaks with single stutter
    SEXP probabilityPeaksS(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                                SEXP DNAcont, SEXP stutterMean, SEXP stutterAdjust, SEXP stutterGradient, 
                                SEXP scale, SEXP degradation, SEXP fragLengths, 
                                SEXP fragNames, SEXP repAdjust, SEXP detectionThresh); 
    //! \brief probability of peaks with single and double stutter
    SEXP probabilityPeaksSD(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                                SEXP DNAcont, SEXP stutterMean, SEXP stutterAdjust, SEXP doubleStutterRate, SEXP stutterGradient, 
                                SEXP scale, SEXP degradation, SEXP fragLengths, 
                                SEXP fragNames, SEXP repAdjust, SEXP detectionThresh);
    //! \brief probability of peaks with single and over stutter
    SEXP probabilityPeaksSO(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                                SEXP DNAcont, SEXP stutterMean, SEXP stutterAdjust, SEXP overStutterRate, SEXP stutterGradient, 
                                SEXP scale, SEXP degradation, SEXP fragLengths, 
                                SEXP fragNames, SEXP repAdjust, SEXP detectionThresh);
    //! \brief probability of peaks with single, double and over stutter
    SEXP probabilityPeaksSDO(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                                SEXP DNAcont, SEXP stutterMean, SEXP stutterAdjust, SEXP doubleStutterRate, SEXP overStutterRate, SEXP stutterGradient, 
                                SEXP scale, SEXP degradation, SEXP fragLengths, 
                                SEXP fragNames, SEXP repAdjust, SEXP detectionThresh);

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

#ifdef __cplusplus
}
#endif
#endif
