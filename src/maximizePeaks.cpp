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






//std::vector<double> peakMeanDose(SEXP genotype,SEXP alleles,SEXP heights,SEXP sizes,SEXP rcont,SEXP stutter,SEXP degradation,SEXP fragLengths, SEXP fragNames)
SEXP peakMeanDose(SEXP genotype, SEXP alleles,SEXP heights,SEXP sizes,SEXP DNAcont,SEXP stutter,SEXP degradation,SEXP fragL, SEXP fragN,SEXP repAdjust)
	{
	int const nGen = length(genotype);
	int const nCSP = length(alleles);
	int const nCont = length(DNAcont);
	int const nFrag = length(fragN);
	SEXP genDuplicate = duplicate(genotype);
	double const * const genotype_ptr     = REAL(genDuplicate);
	double const * const allele_ptr     = REAL(alleles);
	double const * const height_ptr     = REAL(heights);
	double const * const size_ptr     = REAL(sizes);
	double const * const DNAcont_ptr     = REAL(DNAcont);
	double const * const stutter_ptr     = REAL(stutter);
	double const * const deg_ptr     = REAL(degradation);
	double const * const repadjust_ptr     = REAL(repAdjust);
	double const * const fragL_ptr     = REAL(fragL);
	double DNAcontSub, degSub, tmpDose;
	SEXP stutterPos = PROTECT(allocVector(REALSXP, nGen));
	double * stutterPos_ptr     = REAL(stutterPos);
	SEXP allPos = PROTECT(allocVector(REALSXP, 2*nGen));
	double * allPos_ptr     = REAL(allPos);
	SEXP fragNvec = PROTECT(duplicate(fragN));
	double * fragNvec_ptr     = REAL(fragNvec);
	float fragSub;
	std::vector<double> outMu;
	std::vector<float> fragVecN, fragVecL, genotypeVec, stutterPosVec, allPosVec, debug;
	std::vector<float>::iterator itFlt,itFlt2;
	genoStruct tmpMu;
	std::vector<genoStruct> muA, muS, outRes;

    	// convert fragLengths to vector
	for(int i=0; i<nFrag; i++)
		{
		fragVecN.push_back(std::floor(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragL_ptr[i]);
		}

	// Loop over members of genotype
	for(int y=0; y<nGen; y++)
		{
	        // round genotypes
		genotypeVec.push_back(std::floor(genotype_ptr[y]*10.0f)/10.0f);
		// get stutter positions
		stutterPosVec.push_back(std::floor((genotype_ptr[y]-1.0)*10.0f)/10.0f);
		// get all positions, while rounding to 1dp
		allPosVec.push_back(genotypeVec[y]);
		allPosVec.push_back(stutterPosVec[y]);
		}
    	PROTECT(genDuplicate);
    	// sort and unique allPos
		std::sort(allPosVec.begin(),allPosVec.end());
		itFlt = std::unique(allPosVec.begin(),allPosVec.end());
		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; i++)
		{
		// DNAcont and deg
		if(i%2!=0)
			{
			DNAcontSub = DNAcont_ptr[(i-1)/2];
			degSub = deg_ptr[(i-1)/2];
			} else {
			DNAcontSub = DNAcont_ptr[i/2];
			degSub = deg_ptr[i/2];
			}
		// fragLengths
		itFlt2 = std::find(fragVecN.begin(),fragVecN.end(),genotypeVec[i]);
		fragSub = fragVecL[std::distance(fragVecN.begin(),itFlt2)];
		// effective dose
		tmpDose = repadjust_ptr[0]*DNAcontSub*std::pow((1+degSub),fragSub);

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
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes.push_back(tmpMu);
		}
	
	// convert outMu into SEXP object (for debugging etc)
	SEXP outsexp = PROTECT(allocVector(REALSXP, allPosVec.size()));
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
	UNPROTECT(5);
	return outsexp;
	}

struct sortByAlleleNameCSP
    {
    bool operator()(cspStruct const  &a, cspStruct const &b)
        {
        return(a.alleles>b.alleles);
        }
    };

struct sortByAlleleNameGeno
    {
    bool operator()(genoStruct const  &a, genoStruct const &b)
        {
        return(a.genotype>b.genotype);
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

// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
    if (x <= 0 || a <= 0 || b <= 0)
        return 0.0;
    else
        return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}


SEXP probabilityPeaks(SEXP genotypeArray, SEXP alleles, SEXP heights, SEXP sizes, SEXP DNAcont, SEXP stutter, SEXP scale, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP repAdjust)
	{
	cspStruct tmpCSP;
	std::vector<cspStruct> csp, cspModify;
	SEXP gammaMu;
	std::vector<genoStruct> gammaMuVec;
	genoStruct tmpGeno;
	std::vector<genoStruct>::iterator genoIt;
	int nondropoutFlag, unhypFlag;
	bool removeFlag;
	double outDensity, tmpDensity, tmpScale, shape, height, debug;

	genotypeArray = coerceVector(genotypeArray,REALSXP);

	int const nGen = length(DNAcont)*2; // rows
	int const nCombs = length(genotypeArray)/nGen;  // cols
	int const nFrag = length(fragNames);
	int const nCSP = length(alleles);
	int row, column;
	SEXP ALLELES = duplicate(alleles);
	double const * const genotypeArray_ptr     = REAL(genotypeArray);
	double const * const allele_ptr     = REAL(ALLELES);
	double const * const height_ptr     = REAL(heights);
	double const * const size_ptr     = REAL(sizes);
	double const * const scale_ptr     = REAL(scale);
	double const * const fragL_ptr     = REAL(fragLengths);
	double const * const fragN_ptr     = REAL(fragNames);
	std::vector<float> genotypeVec, stutterPosVec, allPosVec, fragVecN, fragVecL;
	std::vector<float>::iterator itFlt;
	std::vector<int> toRemove;
	std::vector<double> debugVec;
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);



    	SEXP tmpVal;
	tmpVal = allocVector(REALSXP,1);
	double *tmp_ptr = REAL(tmpVal);
	
    	// convert fragLengths to vector
	for(int i=0; i<nFrag; i++)
		{
		fragVecN.push_back(std::floor(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragL_ptr[i]);
		}

	// fill csp
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

	SEXP result, input;
	PROTECT(result = allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);

	// sort csp by allele name
    	sort(csp.begin(), csp.end(), sortByAlleleNameCSP());

    	//Rprintf("%i %i\n",nCombs,nGen);

	SEXP genotype = PROTECT(allocVector(REALSXP, nGen));
	double const * const genotype_ptr     = REAL(genotype);

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	for(int x=0; x<nCombs; x++)
		{
		//Rprintf("Before clearing");
		// Clear vectors from last genotype
		genotypeVec.clear();
		stutterPosVec.clear();
		allPosVec.clear();
		cspModify.clear();
		gammaMuVec.clear();
		toRemove.clear();

		// Loop over members of genotype
		for(int y=0; y<nGen; y++)
			{
			// fill genotype vector
			REAL(genotype)[y] = genotypeArray_ptr[(x*nCombs)+y];
		        // round genotypes
			genotypeVec.push_back(std::floor(genotype_ptr[y]*10.0f)/10.0f);
			// get stutter positions
			stutterPosVec.push_back(std::floor((genotype_ptr[y]-1.0)*10.0f)/10.0f);
			// get all positions, while rounding to 1dp
			allPosVec.push_back(genotypeVec[y]);
			allPosVec.push_back(stutterPosVec[y]);
			}
		
		//Rprintf("Before copying");
		// copy csp vector into cspModify vector
		for(int i=0; i<csp.size(); i++)
			{
			cspModify.push_back(csp[i]);
			}

		//SEXP outSexp 
		// sort and unique allPos
		std::sort(allPosVec.begin(),allPosVec.end());
		itFlt = std::unique(allPosVec.begin(),allPosVec.end());
		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));

	
		//Rprintf("Before mean");
    		// get peak gamma means 
		PROTECT(gammaMu = peakMeanDose(genotype, alleles, heights, sizes, DNAcont, stutter, degradation, fragLengths, fragNames, repAdjust));
		double const * const mu_ptr     = REAL(gammaMu);
		for(int i=0; i<length(gammaMu); i++)
			{
			tmpGeno.genotype = allPosVec[i];
			tmpGeno.dose = mu_ptr[i];
			gammaMuVec.push_back(tmpGeno);
			}

		// Add dropout alleles to peak heights, with peak height 0
		for(int i=0; i<allPosVec.size(); i++)
			{
			nondropoutFlag = 0;
			for(int j=0; j<cspModify.size(); j++)
				{
				if(cspModify[j].alleles==allPosVec[i]) nondropoutFlag += 1; 
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
		    tmpScale = scale_ptr[0];
		    shape = gammaMuVec[i].dose/scale_ptr[0];
		    height = cspModify[i].heights;
		    tmpDensity = gammaPdf(height, shape+0.001, 1/(tmpScale+0.001));
//debugVec.push_back(tmpDensity);
		//if(tmpDensity==DBL_MIN) tmpDensity = 0;
		    if(i==0) 
			{
			outDensity = tmpDensity;
			} else {
	    	    	outDensity = outDensity*tmpDensity;
			}
	   	    }
		
		// set out result for this genotype combination
		//out_ptr[x] = outDensity;
		out_ptr[x] = outDensity;
		UNPROTECT(1);
		}

	// debug output vector
	//SEXP debugSEXP = PROTECT(allocVector(REALSXP, debugVec.size()));
	//double const * const debug_ptr     = REAL(debugSEXP);
	//for(int i=0; i<debugVec.size(); i++) REAL(debugSEXP)[i] = debugVec[i];
	//Rprintf("After main loop");
	UNPROTECT(4);
	//return debugSEXP; 
	return result;
	// get into sexp format
//    	SEXP outsexp = allocVector(REALSXP, debugVec.size());
 //   	Rprintf("%i",debugVec.size());
 //   	double * outsexp_ptr     = REAL(outsexp);
//	for(int i=0; i<debugVec.size(); i++)
//		{
//		outsexp_ptr[i] = debugVec[i];
//		}
    	
 //   	return(outsexp);
	}
//*/

















