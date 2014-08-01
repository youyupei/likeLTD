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


struct peakStruct
	{
	double allele;
	double height;
	double size;
	} ;

struct genoStruct
	{
	double genotype;
	double dose;
	} ;

//std::vector<double> peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames)
SEXP peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames)
	{
	int const nGen = INTEGER(GET_DIM(genotype))[0];
	int const nCSP = INTEGER(GET_DIM(alleles))[0];
	int const nCont = INTEGER(GET_DIM(rcont))[0];
	int const nFrag = INTEGER(GET_DIM(fragNames))[0];
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
	double fragSub;
	std::vector<double> stutterPosVec, allPosVec, outMu, fragVec;
	std::vector<double>::iterator itFlt;
	peakStruct tmpPeak;
	std::vector<peakStruct> peaks;
	genoStruct tmpMu;
	std::vector<genoStruct> muA, muS;
	// fill peaks vector
	//for(int i=0; i<nCSP; i++)
	//	{
	//	tmpPeak.allele = allele_ptr[i]; 
	//	tmpPeak.height = height_ptr[i];
	//	tmpPeak.size = size_ptr[i];
	//	peaks.push_back();
	//	}
	// convert fragLengths to vector
	for(int i=0; i<nFrag; i++)
		{
		fragVec.push_back(fragL_ptr[i]);
		}

	for(int i=0; i<nGen; i++)
		{
		// get stutter positions
		//stutterPos_ptr[i] = genotype_ptr[i]-1.0;
		stutterPosVec.push_back(genotype_ptr[i]-1.0);
		// get all positions
		allPosVec.push_back(stutterPosVec[i]);
		allPosVec.push_back(genotype_ptr[i]);
		// get all positions
		//allPos_ptr[i+nGen] = stutterPos_ptr[i];
		//allPos_ptr[i] = genotype_ptr[i];
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
			rcontSub = rcont_ptr[(i+1)/2];
			degSub = deg_ptr[(i+1)/2];
			} else {
			rcontSub = rcont_ptr[i/2];
			degSub = deg_ptr[i/2];
			}
		// fragLengths
		itFlt = std::find(fragVec.begin(),fragVec.end(),genotype_ptr[i]);
		fragSub = fragVec[std::distance(fragVec.begin(),itFlt)];
		// effective dose
		tmpDose = rcontSub*std::pow((1+degSub),fragSub);
		// stutter adjusted effective dose
		tmpMu.genotype = genotype_ptr[i];
		tmpMu.dose = tmpDose * (1-stutter_ptr[0]);
		muA.push_back(tmpMu);
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
	

	SEXP outsexp = allocVector(REALSXP, allPosVec.size());
	double * outsexp_ptr     = REAL(outsexp);
	for(int i=0; i<allPosVec.size(); i++)
		{
		outsexp_ptr[i] = outMu[i];
		}
	//return outMu;
	return outsexp;
	}

//float gammaDensity(float x,float k,float theta)
//	{
//	float out;
//	out = (1/(tgamma(k)*theta^k))*x^(k-1)*exp(-x/theta);
//	return out
//	}
