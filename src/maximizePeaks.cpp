
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


struct sortByAlleleNameCSP
    {
    bool operator()(cspStruct const  &a, cspStruct const &b)
        {
        return(a.alleles<b.alleles);
        }
    };

struct sortByAlleleNameGeno
    {
    bool operator()(genoStruct const  &a, genoStruct const &b)
        {
        return(a.genotype<b.genotype);
        }
    };

inline int myRound(double d)
    {
        return static_cast<int>(d + 0.5);
    }

inline bool shouldBeRemoved( genoStruct g, std::vector<cspStruct> csp) 
	{
	int index;
	bool removeFlag = false;
	for(unsigned int j=0; j<csp.size(); ++j)
	    {
            if(g.genotype==csp[j].alleles) 
                {
                index += 1;
                continue;
                }
	    }
	if(index==0) removeFlag = true;
	return removeFlag;
	}

/* natural log of the gamma function
   gammln as implemented in the
 * first edition of Numerical Recipes in C */
inline double gammln(double xx)
{
    double x, tmp, ser;
    const static double cof[6]={76.18009172947146,    -86.50532032941677,
                                24.01409824083091,    -1.231739572450155,
                                0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.5066282746310005*ser);
}

// natural log of the gamma distribution PDF
// see: http://compbio.mit.edu/spimap/pub/spimap/src/gamma.cpp
inline double gammalog(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}

/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
inline double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

#define KF_GAMMA_EPS 1e-14

// regularized lower incomplete gamma function, by series expansion
// see: http://barricklab.org/hg/breseq/file/9ab26de12a0b/extern/samtools-0.1.15/bcftools/kfunc.c
static double kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}

static double ln_kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return s * log(z) - z - kf_lgamma(s + 1.) + log(sum);
}

inline std::vector<genoStruct> combineDoses(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS)
	{
	genoStruct tmpMu; 
	std::vector<genoStruct> outRes(allPosVec.size(),tmpMu);
	double tmpDose;
	// combine doses at each allelic position
	for(unsigned int i=0; i<allPosVec.size(); ++i)
		{
		tmpDose = 0;
		for(unsigned int j=0; j<muA.size(); ++j)
			{
			if(muA[j].genotype==allPosVec[i]) tmpDose = tmpDose+muA[j].dose;
			if(muS[j].genotype==allPosVec[i]) tmpDose = tmpDose+muS[j].dose;
			}
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes[i] = tmpMu;
		}
	return(outRes);
	}

inline std::vector<genoStruct> peakMeanDose(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec,double stutterTrue, double stutterFalse,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, int nGen, int nCSP, int nCont, int nFrag)
	{
	double tmpDose, fragSub;
	genoStruct tmpMu;
	std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu);
	int matchIndex;
	double diff;

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; ++i)
		{
		// fragLengths
		matchIndex = 0;
		for(unsigned int q=0; q<fragVecN.size(); ++q)
			{
			diff = std::abs(genotypeVec[i]-fragVecN[q]);
			
			if(diff<0.0001) 
			    {
			    matchIndex = q;
                break;
			    }
			}

		fragSub = fragVecL[matchIndex];

		//itDbl = std::find(fragVecN.begin(),fragVecN.end(),std::floor(genotypeVec[i]*10.0f)/10.0f);
		//fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

		// effective dose
		tmpDose = DNAcontVec[i]*std::pow(degVec[i],fragSub);


		// stutter adjusted effective dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = tmpDose * stutterFalse;
		muA[i] = tmpMu;
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = tmpDose * stutterTrue;
		muS[i] = tmpMu;
		}

/*
	// combine doses at each allelic position
	for(int i=0; i<allPosVec.size(); i++)
		{
		tmpDose = 0;
		for(int j=0; j<muA.size(); j++)
			{
			if(muA[j].genotype==allPosVec[i]) tmpDose = tmpDose+muA[j].dose;
			if(muS[j].genotype==allPosVec[i]) tmpDose = tmpDose+muS[j].dose;
			}
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes.push_back(tmpMu);
		}
*/
	std::vector<genoStruct> outRes = combineDoses(allPosVec,muA,muS);
	return outRes;
	}


inline double getDensity(std::vector<genoStruct> gammaMuVec, std::vector<cspStruct> cspModify, double scale, double cdfArg, double pdfArg)
	{
	double outDensity=1, tmpDensity=0;
	unsigned int vecSize = gammaMuVec.size();
    	// get density
	for(unsigned int i=0; i<vecSize; ++i)
		{
		//if(i==gammaMuVec.size()) break;
		if(cspModify[i].heights==0)
			{
			tmpDensity = kf_gammap(gammaMuVec[i].dose/scale,cdfArg);
			} else {
			tmpDensity = exp(gammalog(cspModify[i].heights, gammaMuVec[i].dose/scale, pdfArg));
			}
		//outDensity = outDensity * tmpDensity;
		outDensity *= tmpDensity;
	   	}
	return(outDensity);
	}

inline std::vector<cspStruct> modifyCSP(std::vector<cspStruct> csp,std::vector<float> allPosVec)
	{
	cspStruct tmpCSP;
	tmpCSP.heights=0;
	tmpCSP.sizes=999;
	int nondropoutFlag, initSize=csp.size();
	// Add dropout alleles to peak heights, with peak height 0
	for(unsigned int i=0; i<allPosVec.size(); ++i)
		{
		nondropoutFlag = 0;
		for(unsigned int j=0; j<csp.size(); ++j)
			{
			if((myRound(csp[j].alleles*10.0f)/10.0f)==(myRound(allPosVec[i]*10.0f)/10.0f)) 
				{
				nondropoutFlag = nondropoutFlag + 1; 
				break;	
				}
			}
		if(nondropoutFlag==0)
			{
			tmpCSP.alleles = allPosVec[i];
			csp.push_back(tmpCSP); 
			}
		}
	
	// Sort CSP by allele name
	if(csp.size()!=initSize) sort(csp.begin(), csp.end(), sortByAlleleNameCSP());

	return(csp);
	}


inline double singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutterTrue, double stutterFalse, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, int currentComb, int nGen, int nCSP, int nCont, int nFrag, double cdfArg, double pdfArg)
//std::vector<double> singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutterMean, double stutterAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, double detectionThresh, int currentComb, int nGen, int nCSP, int nCont, int nFrag)
	{
	std::vector<cspStruct> cspModify;
	std::vector<genoStruct> gammaMuVec;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec(nGen,0), stutterPosVec(nGen,0), allPosVec(nGen*2,0);
	std::vector<float>::iterator itFlt;
	//std::vector<double> debug;

	// Loop over members of genotype
	for(int y=0; y<nGen; ++y)
		{
	        // round genotypes
		genotypeVec[y] = myRound(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f;
		// get stutter positions
		stutterPosVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f;
//debug.push_back(((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f);
		// get all positions, while rounding to 1dp
		allPosVec[y] = genotypeVec[y];
		allPosVec[y+nGen] = stutterPosVec[y];
		}
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
		
	//Rprintf("Before mean");
    	// get peak gamma means 
	gammaMuVec = peakMeanDose(genotypeVec, stutterPosVec, allPosVec, csp, DNAcont, stutterTrue, stutterFalse, degradation, fragLengths, fragNames, nGen, nCSP, nCont, nFrag);

/*
	//Rprintf("Before copying");
	// copy csp vector into cspModify vector
	for(int i=0; i<csp.size(); i++)
		{
		cspModify.push_back(csp[i]);
		}

	// Add dropout alleles to peak heights, with peak height 0
	for(int i=0; i<allPosVec.size(); i++)
		{
		nondropoutFlag = 0;
		for(int j=0; j<cspModify.size(); j++)
			{
			if((std::floor(cspModify[j].alleles*10.0f)/10.0f)==(std::floor(allPosVec[i]*10.0f)/10.0f)) nondropoutFlag += 1; 
			}
		if(nondropoutFlag==0)
			{
			tmpCSP.alleles = allPosVec[i];
			tmpCSP.sizes = 999;
			tmpCSP.heights = 0;
			cspModify.push_back(tmpCSP); 
			}
		}
	
	// Sort CSP by allele name
	sort(cspModify.begin(), cspModify.end(), sortByAlleleNameCSP());
*/

// add dropout alleles to the CSP
cspModify = modifyCSP(csp,allPosVec);

/*
// This section removes alleles that require dropin to explain the CSP
// Currently do not allow dropin, so do not need this section
	//Rprintf("Before toRemove");
	// Remove means for alleles that are not being considered		
	for(int i=0; i<gammaMuVec.size(); i++)
		{
		if(shouldBeRemoved(gammaMuVec[i],cspModify))
			{
			toRemove.push_back(i);
			gammaMuVec[i].genotype = 99999;
			}
		}
	//Rprintf("Before toRemove sort");
	sort(gammaMuVec.begin(), gammaMuVec.end(), sortByAlleleNameGeno());
	//Rprintf("Before toRemove popback");
	if(toRemove.size()!=0)
		{
		for(int i=0; i<toRemove.size(); i++)
			{
			gammaMuVec.pop_back();
			}
		}
*/

/*
	//double outDensity=1, tmpDensity;
    	// get density
	for(int i=0; i<gammaMuVec.size(); i++)
		{
		if(i==gammaMuVec.size()) break;
		if(cspModify[i].heights==0)
			{
			tmpDensity = kf_gammap(gammaMuVec[i].dose/scale,detectionThresh/scale);
			} else {
			tmpDensity = exp(gammalog(cspModify[i].heights, gammaMuVec[i].dose/scale, 1/(scale)));
			}
		outDensity = outDensity * tmpDensity;
	   	}
*/

    	// get density
	double outDensity = getDensity(gammaMuVec,cspModify,scale,cdfArg,pdfArg);

//	return(debug);
	return(outDensity);
	}






SEXP probabilityPeaks(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP DNAcont, SEXP stutter, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust, SEXP detectionThresh)
	{
	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	//genotypeArray = coerceVector(genotypeArray,REALSXP);
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCSP = length(alleles);
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
	cspStruct tmpCSP;
	std::vector<cspStruct> csp(nCSP,tmpCSP);
	std::vector<double> outDouble (nCombs);




	// convert genotypeArray to vector
	SEXP GENOTYPEARRAY = PROTECT(duplicate(genotypeArray));
	double const * const genotypeArray_ptr     = REAL(GENOTYPEARRAY);
	std::vector<double> genotypeArrayVec(nCombs*nGen,0);
	for(int i=0; i<nCombs*nGen; ++i)
		{
		genotypeArrayVec[i] = genotypeArray_ptr[i];
		}

	// convert DNAcont to vector
	SEXP DNACONT = PROTECT(duplicate(DNAcont));
	double const * const dnacont_ptr     = REAL(DNACONT);
	std::vector<double> DNAcontVec;
	for(int i=0; i<nCont; ++i)
		{
		DNAcontVec.push_back(dnacont_ptr[i]);
		}	

	// convert stutter to double
	SEXP STUTTER = PROTECT(duplicate(stutter));
	double const * const stutter_ptr     = REAL(STUTTER);
	double stutterDouble = stutter_ptr[0];

	// convert stutter to double
	//SEXP STUTTERADJUST = PROTECT(duplicate(stutterAdjust));
	//double const * const stutterAdjust_ptr     = REAL(STUTTERADJUST);
	//double stutterAdjustDouble = stutterAdjust_ptr[0];

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}

    	// convert fragLengths and fragNames to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	std::vector<double> fragVecN, fragVecL;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(myRound(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		}

    	// convert repAdjust to double
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	double repadjust = repadjust_ptr[0];

    	// convert threshold to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detectionthresh_ptr     = REAL(DETECTIONTHRESH);
	double detectThresh = detectionthresh_ptr[0];

	// convert alleles, heights and sizes to double in csp structure
	SEXP ALLELES = duplicate(alleles);	
	double const * const allele_ptr     = REAL(ALLELES);
	SEXP HEIGHTS = duplicate(heights);
	double const * const height_ptr     = REAL(HEIGHTS);
	SEXP roundedAlleles = PROTECT(allocVector(REALSXP, nCSP));
	double       * const rounded_ptr  = REAL(roundedAlleles);
	for(int i=0; i<nCSP; ++i)
		{
		REAL(roundedAlleles)[i] = myRound(allele_ptr[i]*10.0f)/10.0f;
        	tmpCSP.alleles = rounded_ptr[i];
        	tmpCSP.heights = height_ptr[i];
	       	csp[i] = tmpCSP;   	        	
		}

	// sort csp by allele name
    std::sort(csp.begin(), csp.end(), sortByAlleleNameCSP());

	double stutterTrue = stutterDouble*repadjust;
	double stutterFalse = (1-stutterDouble)*repadjust;

	double cdfArg = detectThresh/scaleDouble;
	double pdfArg = 1/scaleDouble; 

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	# pragma omp parallel for schedule(dynamic)
	for(int x=0; x<nCombs; ++x)
		{
		outDouble[x] = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterTrue, stutterFalse, scaleDouble, degVec, fragVecL, fragVecN, x, nGen, nCSP, nCont, nFrag, cdfArg, pdfArg);
		//outDouble = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterMeanDouble, stutterAdjustDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, detectThresh, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif

/*
	// Make and return debug object
	SEXP debug;
	PROTECT(debug = allocVector(REALSXP, outDouble.size()));
  	double       * const debug_ptr  = REAL(debug);
	for(int i=0; i<outDouble.size(); i++)
		{
		debug_ptr[i] = outDouble[i];
		}
	UNPROTECT(12);
	return debug;
*/


	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
	UNPROTECT(11);
	return result;

	}







