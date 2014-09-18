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

bool shouldBeRemoved( genoStruct g, std::vector<cspStruct> csp) 
	{
	int index;
	bool removeFlag = false;
	for(int j=0; j<csp.size(); j++)
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

// Lanczos approximation to the gamma function. 
// 
// found on http://www.rskey.org/gamma.htm   
//
double gamm(double x) 
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*M_PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}

/* natural log of the gamma function
   gammln as implemented in the
 * first edition of Numerical Recipes in C */
double gammln(double xx)
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

// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}

// natural log of the gamma distribution PDF
double gammalog(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return -x * b + (a - 1.0) * log(x) + a * log(b) - gammln(a);
}





std::vector<double> peakMeanDose(std::vector<float> genotypeVec, std::vector<float> stutterPosVec, std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec,double stutter,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN,double repAdjust, int nGen, int nCSP, int nCont, int nFrag)
	{
	double DNAcontSub, degSub, tmpDose;
	float fragSub;
	std::vector<double> outMu;
	std::vector<float>::iterator itFlt,itFlt2;
	std::vector<double>::iterator itDbl; 
	genoStruct tmpMu;
	std::vector<genoStruct> muA, muS, outRes;
	std::vector<double> debug;

for(int i=0; i<fragVecN.size(); i++) debug.push_back(fragVecN[i]);

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; i++)
		{
		// DNAcont and deg
		if(i%2!=0)
			{
			DNAcontSub = DNAcontVec[(i-1)/2];
			degSub = degVec[(i-1)/2];
			} else {
			DNAcontSub = DNAcontVec[i/2];
			degSub = degVec[i/2];
			}
		// fragLengths
		itDbl = std::find(fragVecN.begin(),fragVecN.end(),genotypeVec[i]);
		fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

		// effective dose
		tmpDose = repAdjust*DNAcontSub*std::pow((1+degSub),fragSub);

		// stutter adjusted effective dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = tmpDose * (1-stutter);
		muA.push_back(tmpMu);
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = tmpDose * stutter;
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
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes.push_back(tmpMu);
		}
	//return(debug);
	return outMu;
	}


double singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutter, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, int currentComb, int nGen, int nCSP, int nCont, int nFrag)
//std::vector<double> singleGenotype(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double stutter, double scale, std::vector<double> degradation, std::vector<double> fragLengths, std::vector<double> fragNames, double repAdjust, int currentComb, int nGen, int nCSP, int nCont, int nFrag)
	{
	cspStruct tmpCSP;
	std::vector<cspStruct> cspModify;
	genoStruct tmpGeno;
	std::vector<genoStruct> gammaMuVec;
	double outDensity, tmpDensity, tmpScale, shape, height;
	int nondropoutFlag;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec, stutterPosVec, allPosVec;
	std::vector<float>::iterator itFlt;
	std::vector<double> gammaMu;
	std::vector<double> debug;

	// Loop over members of genotype
	for(int y=0; y<nGen; y++)
		{
	        // round genotypes
		genotypeVec.push_back(std::floor(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f);
		// get stutter positions
		stutterPosVec.push_back(std::floor((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f);
		// get all positions, while rounding to 1dp
		allPosVec.push_back(genotypeVec[y]);
		allPosVec.push_back(stutterPosVec[y]);
		}
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
		
	//Rprintf("Before copying");
	// copy csp vector into cspModify vector
	for(int i=0; i<csp.size(); i++)
		{
		cspModify.push_back(csp[i]);
		}
	
	//Rprintf("Before mean");
    	// get peak gamma means 
	gammaMu = peakMeanDose(genotypeVec, stutterPosVec, allPosVec, csp, DNAcont, stutter, degradation, fragLengths, fragNames, repAdjust, nGen, nCSP, nCont, nFrag);

	for(int i=0; i<gammaMu.size(); i++)
		{
		tmpGeno.genotype = allPosVec[i];
		tmpGeno.dose = gammaMu[i];
		gammaMuVec.push_back(tmpGeno);
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



	//Rprintf("Before density");
    	// get density
	for(int i=0; i<gammaMuVec.size(); i++)
		{
		if(i==gammaMuVec.size()) break;
		tmpScale = scale;
		shape = gammaMuVec[i].dose/tmpScale;
		height = cspModify[i].heights;
		tmpDensity = exp(gammalog(height+0.1, shape+0.1, 1/(tmpScale+0.1)));


		if(i==0) 
			{
			outDensity = tmpDensity;
			} else {
	    	outDensity = outDensity * tmpDensity;
			}
	   	}

	debug.push_back(outDensity);
	return(outDensity);
	//return debug;
	}






SEXP probabilityPeaks(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP sizes, SEXP DNAcont, SEXP stutter, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust)
	{
	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	cspStruct tmpCSP;
	std::vector<cspStruct> csp;
	//genotypeArray = coerceVector(genotypeArray,REALSXP);
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCSP = length(alleles);
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
	std::vector<double> outDouble (nCombs);




	// convert genotypeArray to vector
	SEXP GENOTYPEARRAY = PROTECT(duplicate(genotypeArray));
	double const * const genotypeArray_ptr     = REAL(GENOTYPEARRAY);
	std::vector<double> genotypeArrayVec;
	for(int i=0; i<nCombs*nGen; i++)
		{
		genotypeArrayVec.push_back(genotypeArray_ptr[i]);
		}

	// convert DNAcont to vector
	SEXP DNACONT = PROTECT(duplicate(DNAcont));
	double const * const dnacont_ptr     = REAL(DNACONT);
	std::vector<double> DNAcontVec;
	for(int i=0; i<nCont; i++)
		{
		DNAcontVec.push_back(dnacont_ptr[i]);
		}	

	// convert stutter to double
	SEXP STUTTER = PROTECT(duplicate(stutter));
	double const * const stutter_ptr     = REAL(STUTTER);
	double stutterDouble = stutter_ptr[0];

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; i++)
		{
		degVec.push_back(deg_ptr[i]);
		}

    	// convert fragLengths and fragNames to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	std::vector<double> fragVecN, fragVecL;
	for(int i=0; i<nFrag; i++)
		{
		fragVecN.push_back(std::floor(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		}

    	// convert repAdjust to double
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	double repadjust = repadjust_ptr[0];

	// convert alleles, heights and sizes to double in csp structure
	SEXP ALLELES = duplicate(alleles);	
	double const * const allele_ptr     = REAL(ALLELES);
	SEXP HEIGHTS = duplicate(heights);
	double const * const height_ptr     = REAL(HEIGHTS);
	SEXP SIZES = duplicate(sizes);
	double const * const size_ptr     = REAL(SIZES);
	SEXP roundedAlleles = PROTECT(allocVector(REALSXP, nCSP));
	double       * const rounded_ptr  = REAL(roundedAlleles);
	for(int i=0; i<nCSP; i++)
		{
		REAL(roundedAlleles)[i] = std::floor(allele_ptr[i]*10.0f)/10.0f;
        	tmpCSP.alleles = rounded_ptr[i];
        	tmpCSP.heights = height_ptr[i];
        	tmpCSP.sizes = size_ptr[i];
	       	csp.push_back(tmpCSP);   	        	
		}

	// sort csp by allele name
    	std::sort(csp.begin(), csp.end(), sortByAlleleNameCSP());

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	int x=0;
	# pragma omp parallel for //private(genotypeArrayVec, csp, DNAcontVec, stutterDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, x)
	for(x=0; x<nCombs; x++)
		{
		outDouble[x] = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif


	// Make and return output object
	SEXP result, input;
	PROTECT(result = allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; i++)
		{
//Rprintf("%f \n",outDouble[i]);
		out_ptr[i] = outDouble[i];
		}
	UNPROTECT(10);
	return result;

	}








