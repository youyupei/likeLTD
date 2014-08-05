#include "config.h"
#include "openmp.h"
#include "maximizePeaks.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <R_ext/Error.h>


#ifdef _OPENMP
#  include <omp.h>
#  ifdef OPENMP_STACK
#    define CSTACK_DEFNS 7
//#    include <Rinterface.h>
#  endif
#endif


struct genoStruct
	{
	double genotype;
	double dose;
	} ;

//std::vector<double> peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames)
SEXP peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames)
	{
	int const nGen = length(genotype);
	int const nCSP = length(alleles);
	int const nCont = length(rcont);
	int const nFrag = length(fragNames);
	double const * const genotype_ptr     = REAL(genotype);
	double const * const allele_ptr     = REAL(alleles);
	double const * const height_ptr     = REAL(heights);
	double const * const size_ptr     = REAL(sizes);
	double const * const rcont_ptr     = REAL(rcont);
	double const * const stutter_ptr     = REAL(stutter);
	double const * const deg_ptr     = REAL(degradation);
	double const * const fragL_ptr     = REAL(fragLengths);
	double const * const fragN_ptr     = REAL(fragNames);
	double rcontSub, degSub, tmpDose;
	//SEXP stutterPos = allocVector(REALSXP, nGen);
	//double * stutterPos_ptr     = REAL(stutterPos);
	//SEXP allPos = allocVector(REALSXP, 2*nGen);
	//double * allPos_ptr     = REAL(allPos);
	float fragSub;
	std::vector<double> outMu;
	std::vector<float> genotypeVec, stutterPosVec, allPosVec, fragVecN, fragVecL, debug;
	std::vector<float>::iterator itFlt,itFlt2;
	genoStruct tmpMu;
	std::vector<genoStruct> muA, muS;
	// convert fragLengths to vector
	for(int i=0; i<nFrag; i++)
		{
		fragVecN.push_back(std::floor(fragN_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragL_ptr[i]);
		}

	for(int i=0; i<nGen; i++)
		{
		// round genotypes
		genotypeVec.push_back(std::floor(genotype_ptr[i]*10.0f)/10.0f);
		// get stutter positions
		stutterPosVec.push_back(std::floor((genotypeVec[i]-1.0)*10.0f)/10.0f);
		// get all positions, while rounding to 1dp
		allPosVec.push_back(stutterPosVec[i]);
		allPosVec.push_back(genotypeVec[i]);
		}
	//SEXP outSexp 
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
	// effective dose for stutter and allelic
	for(int i=0; i<nGen; i++)
		{
		// rcont and deg
		if(i%2!=0)
			{
			rcontSub = rcont_ptr[(i-1)/2];
			degSub = deg_ptr[(i-1)/2];
			} else {
			rcontSub = rcont_ptr[i/2];
			degSub = deg_ptr[i/2];
			}
		// fragLengths
		itFlt2 = std::find(fragVecN.begin(),fragVecN.end(),genotypeVec[i]);
		fragSub = fragVecL[std::distance(fragVecN.begin(),itFlt2)];
		// effective dose
		tmpDose = rcontSub*std::pow((1+degSub),fragSub);

		// stutter adjusted effective dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = tmpDose * (1-stutter_ptr[0]);
		muA.push_back(tmpMu);
		debug.push_back(tmpMu.dose);
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = tmpDose * stutter_ptr[0];
		muS.push_back(tmpMu);

		}
	// combine doses at each allelic position
	for(int i=0; i<allPosVec.size(); i++)
		{
		tmpDose = 0;
		for(int j=0; j<muA.size(); j++)
			{
			if(muA[j].genotype==allPosVec[i]) tmpDose = tmpDose+muA[j].dose;
			if(muS[j].genotype==allPosVec[i]) tmpDose = tmpDose+muS[j].dose;
			}
		outMu.push_back(tmpDose);
		}
	
	// convert outMu into SEXP object (for debugging etc)
	SEXP outsexp = allocVector(REALSXP, allPosVec.size());
	double * outsexp_ptr     = REAL(outsexp);
	for(int i=0; i<allPosVec.size(); i++)
		{
		outsexp_ptr[i] = outMu[i];
		}
//	for(int i=0; i<debug.size(); i++)
//		{
//		outsexp_ptr[i] = debug[i];
//		}
	//return outMu;
	return outsexp;
	}

//float gammaDensity(float x,float k,float theta)
//	{
//	float out;
//	out = (1/(tgamma(k)*theta^k))*x^(k-1)*exp(-x/theta);
//	return out
//	}
