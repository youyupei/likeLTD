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
	        float genotype;
	        float dose;
	        } ;

	//! \brief structure to contain all data associated with a CSP peak
	struct cspStruct
		{
		float alleles;
		float heights;
		float sizes;
		} ;
	
	//! \brief peak mean dose
        //SEXP peakMeanDose(SEXP genotype, SEXP alleles,SEXP heights,SEXP sizes,SEXP DNAcont,SEXP stutter,SEXP degradation,SEXP fragL, SEXP fragN,SEXP repAdjust);
	std::vector<double> peakMeanDose(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec,double stutter,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN,double repAdjust, int nGen, int nCSP, int nCont, int nFrag);

    	//! \brief probability of peaks
        SEXP probabilityPeaks(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP sizes, SEXP DNAcont, SEXP stutter, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust);

	//! \brief probability of a single genotype
	double singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutter, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, int currentComb, int nGen, int nCSP, int nCont, int nFrag);
	//std::vector<double> singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutter, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, int currentComb, int nGen, int nCSP, int nCont, int nFrag);
	
	//! \brief gamma function
	double gamm(double x);

	//! \brief log gamma function
	double gammln(double xx);

	//! \brief pdf of gamma distribution
	double gammaPdf(double x, double a, double b);

	//! \brief log pdf of gamma distribution
	double gammalog(double x, double a, double b);

    //! \brief function to know which of gammaMu to not include
	bool shouldBeRemoved( genoStruct g, std::vector<cspStruct> csp);
#ifdef __cplusplus
}
#endif
#endif
