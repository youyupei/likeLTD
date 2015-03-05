
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

// combine doses, taking into account Single stutter
inline std::vector<genoStruct> combineDosesS(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS)
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

// combine doses taking into account Single and Double stutter 
inline std::vector<genoStruct> combineDosesSD(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS
                                            ,std::vector<genoStruct> muSd)
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
			if(muSd[j].genotype==allPosVec[i]) tmpDose = tmpDose+muSd[j].dose;
			}
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes[i] = tmpMu;
		}
	return(outRes);
	}

// combine doses taking into account Single and Over stuttter
inline std::vector<genoStruct> combineDosesSO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS
                                            ,std::vector<genoStruct> muSo)
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
			if(muSo[j].genotype==allPosVec[i]) tmpDose = tmpDose+muSo[j].dose;
			}
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes[i] = tmpMu;
		}
	return(outRes);
	}

// combine doses taking into account Single, Double and Over stutter
inline std::vector<genoStruct> combineDosesSDO(std::vector<float> allPosVec,std::vector<genoStruct> muA,std::vector<genoStruct> muS
                                            ,std::vector<genoStruct> muSd,std::vector<genoStruct> muSo)
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
			if(muSd[j].genotype==allPosVec[i]) tmpDose = tmpDose+muSd[j].dose;
			if(muSo[j].genotype==allPosVec[i]) tmpDose = tmpDose+muSo[j].dose;
			}
		tmpMu.genotype = allPosVec[i];
		tmpMu.dose = tmpDose;
		outRes[i] = tmpMu;
		}
	return(outRes);
	}

// get expected peak heights taking into account Single stutter
inline std::vector<genoStruct> peakMeanDoseS(std::vector<float> genotypeVec, std::vector<float> stutterPosVec,
                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec, double gradientS, 
                                            double meanS, 
                                            double repAdjust,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
                                            std::vector<double> stutterIndex, int nGen, int nCSP, int nCont, int nFrag)
	{
	double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
	genoStruct tmpMu;
	std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu);
	int matchIndex;
	double diff;

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; ++i)
		{
		// Do not compute dose again if homozygote
		if((genotypeVec[i]==genotypeVec[i-1])&&(i%2!=0)&&(i>0))
			{
			
			} else {
			// get fragLengths
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
			stuttIndSub = stutterIndex[matchIndex];

		//itDbl = std::find(fragVecN.begin(),fragVecN.end(),std::floor(genotypeVec[i]*10.0f)/10.0f);
		//fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

			// compute effective dose
			tmpDose = DNAcontVec[i]*std::pow(degVec[i],fragSub)*repAdjust;
			//double stutterRate = meanS*(1+gradientS*stuttIndSub);
			double stutterRate = meanS+(gradientS*stuttIndSub);
			//stutterDose = tmpDose * stutter;
			stutterDose = tmpDose * stutterRate;
			//nonstutterDose = tmpDose * (1-stutter);
			nonstutterDose = tmpDose * (1-stutterRate);

			}

		// stutter adjusted effective dose
		// non-stutter dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = nonstutterDose;
		muA[i] = tmpMu;
		// stutter dose
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = stutterDose;
		muS[i] = tmpMu;
		}
	std::vector<genoStruct> outRes = combineDosesS(allPosVec,muA,muS);
	return outRes;
	}

// get expected peak heights taking into account Single and Double stutter
inline std::vector<genoStruct> peakMeanDoseSD(std::vector<float> genotypeVec, std::vector<float> stutterPosVec,
                                            std::vector<float> doubleStutterVec,
                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec, double gradientS, 
                                            double meanD,
                                            double meanS, 
                                            double repAdjust,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
                                            std::vector<double> stutterIndex, int nGen, int nCSP, int nCont, int nFrag)
	{
	double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
	genoStruct tmpMu;
	std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu);
	int matchIndex;
	double diff;

	double doubleStutterDose;

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; ++i)
		{
		// Do not compute dose again if homozygote
		if((genotypeVec[i]==genotypeVec[i-1])&&(i%2!=0)&&(i>0))
			{
			
			} else {
			// get fragLengths
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
			stuttIndSub = stutterIndex[matchIndex];

		//itDbl = std::find(fragVecN.begin(),fragVecN.end(),std::floor(genotypeVec[i]*10.0f)/10.0f);
		//fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

			// compute effective dose
			tmpDose = DNAcontVec[i]*std::pow(degVec[i],fragSub)*repAdjust;
			//double stutterRate = meanS*(1+gradientS*stuttIndSub);
			double stutterRate = meanS+(gradientS*stuttIndSub);
			//stutterDose = tmpDose * stutter;
			stutterDose = tmpDose * stutterRate;
			doubleStutterDose = tmpDose * meanD;
			//stutterDose = tmpDose * stutterRate;
			//nonstutterDose = tmpDose * (1-stutter);
			//nonstutterDose = tmpDose * (1-(stutter+doubleStutterRate));
			nonstutterDose = tmpDose * (1-(stutterRate+meanD));
			//nonstutterDose = tmpDose * (1-stutterRate);

			}

		// stutter adjusted effective dose
		// non-stutter dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = nonstutterDose;
		muA[i] = tmpMu;
		// stutter dose
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = stutterDose;
		muS[i] = tmpMu;
		// double stutter dose
		tmpMu.genotype = doubleStutterVec[i];
		tmpMu.dose = doubleStutterDose;
		muSd[i] = tmpMu;
		}
	std::vector<genoStruct> outRes = combineDosesSD(allPosVec,muA,muS,muSd);
	return outRes;
	}

// get expected peak heights taking into account Single and Over stutter
inline std::vector<genoStruct> peakMeanDoseSO(std::vector<float> genotypeVec, std::vector<float> stutterPosVec,
                                             std::vector<float> overStutterVec,
                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec, double gradientS, 
                                             double meanO,
                                            double meanS, 
                                            double repAdjust,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
                                            std::vector<double> stutterIndex, int nGen, int nCSP, int nCont, int nFrag)
	{
	double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
	genoStruct tmpMu;
	std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
	int matchIndex;
	double diff;

	double overStutterDose;

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; ++i)
		{
		// Do not compute dose again if homozygote
		if((genotypeVec[i]==genotypeVec[i-1])&&(i%2!=0)&&(i>0))
			{
			
			} else {
			// get fragLengths
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
			stuttIndSub = stutterIndex[matchIndex];

		//itDbl = std::find(fragVecN.begin(),fragVecN.end(),std::floor(genotypeVec[i]*10.0f)/10.0f);
		//fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

			// compute effective dose
			tmpDose = DNAcontVec[i]*std::pow(degVec[i],fragSub)*repAdjust;
			//double stutterRate = meanS*(1+gradientS*stuttIndSub);
			double stutterRate = meanS+(gradientS*stuttIndSub);
			//stutterDose = tmpDose * stutter;
			stutterDose = tmpDose * stutterRate;
			overStutterDose = tmpDose * meanO;
			//stutterDose = tmpDose * stutterRate;
			//nonstutterDose = tmpDose * (1-stutter);
			//nonstutterDose = tmpDose * (1-(stutter+doubleStutterRate));
			nonstutterDose = tmpDose * (1-(stutterRate+meanO));
			//nonstutterDose = tmpDose * (1-stutterRate);

			}

		// stutter adjusted effective dose
		// non-stutter dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = nonstutterDose;
		muA[i] = tmpMu;
		// stutter dose
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = stutterDose;
		muS[i] = tmpMu;
		// over stutter dose
		tmpMu.genotype = overStutterVec[i];
		tmpMu.dose = overStutterDose;
		muSo[i] = tmpMu;
		}
	std::vector<genoStruct> outRes = combineDosesSO(allPosVec,muA,muS,muSo);
	return outRes;
	}

// get expected peak heights taking into account Single, Double and Over stutter
inline std::vector<genoStruct> peakMeanDoseSDO(std::vector<float> genotypeVec, std::vector<float> stutterPosVec,
                                            std::vector<float> doubleStutterVec, std::vector<float> overStutterVec,
                                            std::vector<float> allPosVec, std::vector<cspStruct> csp, std::vector<double> DNAcontVec, double gradientS, 
                                            double meanD, double meanO,
                                            double meanS, 
                                            double repAdjust,std::vector<double> degVec,std::vector<double> fragVecL, std::vector<double> fragVecN, 
                                            std::vector<double> stutterIndex, int nGen, int nCSP, int nCont, int nFrag)
	{
	double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
	genoStruct tmpMu;
	std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
	int matchIndex;
	double diff;

	double doubleStutterDose, overStutterDose;

	// effective dose for stutter and allelic
	for(int i=0; i<nGen; ++i)
		{
		// Do not compute dose again if homozygote
		if((genotypeVec[i]==genotypeVec[i-1])&&(i%2!=0)&&(i>0))
			{
			
			} else {
			// get fragLengths
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
			stuttIndSub = stutterIndex[matchIndex];

		//itDbl = std::find(fragVecN.begin(),fragVecN.end(),std::floor(genotypeVec[i]*10.0f)/10.0f);
		//fragSub = fragVecL[std::distance(fragVecN.begin(),itDbl)];

			// compute effective dose
			tmpDose = DNAcontVec[i]*std::pow(degVec[i],fragSub)*repAdjust;
			//double stutterRate = meanS*(1+gradientS*stuttIndSub);
			double stutterRate = meanS+(gradientS*stuttIndSub);
			//stutterDose = tmpDose * stutter;
			stutterDose = tmpDose * stutterRate;
			doubleStutterDose = tmpDose * meanD;
			overStutterDose = tmpDose * meanO;
			//stutterDose = tmpDose * stutterRate;
			//nonstutterDose = tmpDose * (1-stutter);
			//nonstutterDose = tmpDose * (1-(stutter+doubleStutterRate));
			nonstutterDose = tmpDose * (1-(stutterRate+meanD+meanO));
			//nonstutterDose = tmpDose * (1-stutterRate);

			}

		// stutter adjusted effective dose
		// non-stutter dose
		tmpMu.genotype = genotypeVec[i];
		tmpMu.dose = nonstutterDose;
		muA[i] = tmpMu;
		// stutter dose
		tmpMu.genotype = stutterPosVec[i];
		tmpMu.dose = stutterDose;
		muS[i] = tmpMu;
		// double stutter dose
		tmpMu.genotype = doubleStutterVec[i];
		tmpMu.dose = doubleStutterDose;
		muSd[i] = tmpMu;
		// over stutter dose
		tmpMu.genotype = overStutterVec[i];
		tmpMu.dose = overStutterDose;
		muSo[i] = tmpMu;
		}
	std::vector<genoStruct> outRes = combineDosesSDO(allPosVec,muA,muS,muSd,muSo);
	return outRes;
	}




inline double getDensity(std::vector<genoStruct> gammaMuVec, std::vector<cspStruct> cspModify, double scale, double cdfArg, double pdfArg)
	{
	double outDensity=0, tmpDensity=0;
	//double outDensity=1, tmpDensity=0;
	unsigned int vecSize = gammaMuVec.size();
    	// get density
	for(unsigned int i=0; i<vecSize; ++i)
		{
		//if(i==gammaMuVec.size()) break;
		if(cspModify[i].heights==0)
			{
			tmpDensity = ln_kf_gammap(gammaMuVec[i].dose/scale,cdfArg);
			//tmpDensity = kf_gammap(gammaMuVec[i].dose/scale,cdfArg);
			} else {
			tmpDensity = gammalog(cspModify[i].heights, gammaMuVec[i].dose/scale, pdfArg);
			//tmpDensity = exp(gammalog(cspModify[i].heights, gammaMuVec[i].dose/scale, pdfArg));
			}
		outDensity += tmpDensity;
		//outDensity *= tmpDensity;
	   	}
	return(exp(outDensity));
	//return(outDensity);
	}

inline std::vector<cspStruct> modifyCSP(std::vector<cspStruct> csp,std::vector<float> allPosVec)
	{
	cspStruct tmpCSP;
	tmpCSP.heights=0;
	tmpCSP.sizes=999;
	int nondropoutFlag, initSize=csp.size();
	//int nondropoutFlag;
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
	//sort(csp.begin(), csp.end(), sortByAlleleNameCSP());

	return(csp);
	}


// Get probability of single genotype combination, taking into account Single stutter
inline double singleGenotypeS(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double gradientS, 
                            double meanS, 
                            double repAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
                            std::vector<double> fragNames, std::vector<double> stutterIndex, int currentComb, int nGen, int nCSP, int nCont, int nFrag, 
                            double cdfArg, double pdfArg)
	{
	std::vector<cspStruct> cspModify;
	std::vector<genoStruct> gammaMuVec;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec(nGen,0), stutterPosVec(nGen,0), allPosVec(nGen*2,0);;
	std::vector<float>::iterator itFlt;
	std::vector<double>::iterator itDbl;
	double gen;

	// Loop over members of genotype
	for(int y=0; y<nGen; ++y)
		{
	    // round genotypes
		genotypeVec[y] = myRound(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f;
		// get stutter positions
		stutterPosVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f;
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
	gammaMuVec = peakMeanDoseS(genotypeVec, stutterPosVec, 
	                        allPosVec, csp, DNAcont, gradientS, 
	                        meanS, 
	                        repAdjust, degradation, fragLengths, 
	                        fragNames, stutterIndex, nGen, nCSP, nCont, nFrag);

    // add dropout alleles to the CSP
    cspModify = modifyCSP(csp,allPosVec);


    // get density
	double outDensity = getDensity(gammaMuVec,cspModify,scale,cdfArg,pdfArg);

    //	return(debug);
	return(outDensity);
	}

// Get probability of single genotype combination, taking into account Single and Double stutter
inline double singleGenotypeSD(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double gradientS, 
                            double meanD,
                            double meanS, 
                            double repAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
                            std::vector<double> fragNames, std::vector<double> stutterIndex, int currentComb, int nGen, int nCSP, int nCont, int nFrag, 
                            double cdfArg, double pdfArg)
	{
	std::vector<cspStruct> cspModify;
	std::vector<genoStruct> gammaMuVec;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), allPosVec(nGen*3,0);
	std::vector<float>::iterator itFlt;
	std::vector<double>::iterator itDbl;
	double gen;


	// Loop over members of genotype
	for(int y=0; y<nGen; ++y)
		{
	    // round genotypes
		genotypeVec[y] = myRound(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f;
		// get stutter positions
		stutterPosVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f;
	    // double stutter positions
        doubleStutterVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-2.0)*10.0f)/10.0f;
		// get all positions, while rounding to 1dp
		allPosVec[y] = genotypeVec[y];
		allPosVec[y+nGen] = stutterPosVec[y];
		// double stutter
		allPosVec[y+(nGen*2)] = doubleStutterVec[y];
		}
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
		
	//Rprintf("Before mean");
    // get peak gamma means 
	gammaMuVec = peakMeanDoseSD(genotypeVec, stutterPosVec, 
	                                        doubleStutterVec,
	                                        allPosVec, csp, DNAcont, gradientS, 
	                                        meanD,
	                                        meanS, 
	                                        repAdjust, degradation, fragLengths, 
	                                        fragNames, stutterIndex, nGen, nCSP, nCont, nFrag);


    // add dropout alleles to the CSP
    cspModify = modifyCSP(csp,allPosVec);


    // get density
	double outDensity = getDensity(gammaMuVec,cspModify,scale,cdfArg,pdfArg);

    //	return(debug);
	return(outDensity);
	}

// Get probability of single genotype combination, taking into account Single and Over stutter
inline double singleGenotypeSO(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double gradientS, 
                            double meanO,
                            double meanS, 
                            double repAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
                            std::vector<double> fragNames, std::vector<double> stutterIndex, int currentComb, int nGen, int nCSP, int nCont, int nFrag, 
                            double cdfArg, double pdfArg)
	{
	std::vector<cspStruct> cspModify;
	std::vector<genoStruct> gammaMuVec;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec(nGen,0), stutterPosVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*3,0);
	std::vector<float>::iterator itFlt;
	std::vector<double>::iterator itDbl;
	double gen;


	// Loop over members of genotype
	for(int y=0; y<nGen; ++y)
		{
	    // round genotypes
		genotypeVec[y] = myRound(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f;
		// get stutter positions
		stutterPosVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f;
	    // over stutter positions
        overStutterVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]+1.0)*10.0f)/10.0f;
		// get all positions, while rounding to 1dp
		allPosVec[y] = genotypeVec[y];
		allPosVec[y+nGen] = stutterPosVec[y];
		// over stutter
		allPosVec[y+(nGen*2)] = overStutterVec[y];
		}
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
		
	//Rprintf("Before mean");
    // get peak gamma means 
	gammaMuVec = peakMeanDoseSO(genotypeVec, stutterPosVec, 
	                                        overStutterVec,
	                                        allPosVec, csp, DNAcont, gradientS, 
	                                        meanO,
	                                        meanS, 
	                                        repAdjust, degradation, fragLengths, 
	                                        fragNames, stutterIndex, nGen, nCSP, nCont, nFrag);


    // add dropout alleles to the CSP
    cspModify = modifyCSP(csp,allPosVec);


    // get density
	double outDensity = getDensity(gammaMuVec,cspModify,scale,cdfArg,pdfArg);

    //	return(debug);
	return(outDensity);
	}

// Get probability of single genotype combination, taking into account Single, Double and Over stutter
inline double singleGenotypeSDO(std::vector<double> genotypeArray, std::vector<cspStruct> csp, std::vector<double> DNAcont, double gradientS, 
                            double meanD, double meanO,
                            double meanS, 
                            double repAdjust, double scale, std::vector<double> degradation, std::vector<double> fragLengths, 
                            std::vector<double> fragNames, std::vector<double> stutterIndex, int currentComb, int nGen, int nCSP, int nCont, int nFrag, 
                            double cdfArg, double pdfArg)
	{
	std::vector<cspStruct> cspModify;
	std::vector<genoStruct> gammaMuVec;
	std::vector<int> toRemove;
	std::vector<float> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*4,0);
	std::vector<float>::iterator itFlt;
	std::vector<double>::iterator itDbl;
	double gen;


	// Loop over members of genotype
	for(int y=0; y<nGen; ++y)
		{
	    // round genotypes
		genotypeVec[y] = myRound(genotypeArray[(currentComb*nGen)+y]*10.0f)/10.0f;
		// get stutter positions
		stutterPosVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-1.0)*10.0f)/10.0f;
	    // double stutter positions
        doubleStutterVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]-2.0)*10.0f)/10.0f;
	    // over stutter positions
        overStutterVec[y] = myRound((genotypeArray[(currentComb*nGen)+y]+1.0)*10.0f)/10.0f;
		// get all positions, while rounding to 1dp
		allPosVec[y] = genotypeVec[y];
		allPosVec[y+nGen] = stutterPosVec[y];
		// double stutter
		allPosVec[y+(nGen*2)] = doubleStutterVec[y];
		// over stutter
		allPosVec[y+(nGen*3)] = overStutterVec[y];
		}
	// sort and unique allPos
	std::sort(allPosVec.begin(),allPosVec.end());
	itFlt = std::unique(allPosVec.begin(),allPosVec.end());
	allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
		
	//Rprintf("Before mean");
    // get peak gamma means 
	gammaMuVec = peakMeanDoseSDO(genotypeVec, stutterPosVec, 
	                                        doubleStutterVec,overStutterVec,
	                                        allPosVec, csp, DNAcont, gradientS, 
	                                        meanD,meanO,
	                                        meanS, 
	                                        repAdjust, degradation, fragLengths, 
	                                        fragNames, stutterIndex, nGen, nCSP, nCont, nFrag);


    // add dropout alleles to the CSP
    cspModify = modifyCSP(csp,allPosVec);


    // get density
	double outDensity = getDensity(gammaMuVec,cspModify,scale,cdfArg,pdfArg);

    //	return(debug);
	return(outDensity);
	}



// get probability of all genotype combinations, taking into account Single and Double stutter
SEXP probabilityPeaksS(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                    SEXP DNAcont, SEXP gradientS,
                    SEXP meanS, 
                    SEXP scale, SEXP degradation, SEXP fragLengths, 
                    SEXP fragNames, SEXP stutterIndex, SEXP repAdjust, SEXP detectionThresh)
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
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter to double
	SEXP MEANS = PROTECT(duplicate(meanS));
	double const * const meanS_ptr     = REAL(MEANS);
	double means = meanS_ptr[0];

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
	SEXP stutterindex = PROTECT(duplicate(stutterIndex));
	double * stutterindex_ptr     = REAL(stutterindex);
	std::vector<double> fragVecN, fragVecL,indexStutter;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(myRound(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		indexStutter.push_back(stutterindex_ptr[i]);
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

	double cdfArg = detectThresh/scaleDouble;
	double pdfArg = 1/scaleDouble; 

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	# pragma omp parallel for schedule(dynamic)
	for(int x=0; x<nCombs; ++x)
		{
		outDouble[x] = singleGenotypeS(genotypeArrayVec, csp, DNAcontVec, gradients, 
		means, 
		repadjust, scaleDouble, degVec, fragVecL, fragVecN, indexStutter, x, nGen, nCSP, nCont, nFrag, cdfArg, pdfArg);
		//outDouble = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterMeanDouble, stutterAdjustDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, detectThresh, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif

	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}

    UNPROTECT(13);
	return result;

	}

// get probability of all genotype combinations, taking into account Single and Double stutter
SEXP probabilityPeaksSD(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                    SEXP DNAcont, SEXP gradientS, SEXP meanD,
                    SEXP meanS, 
                    SEXP scale, SEXP degradation, SEXP fragLengths, 
                    SEXP fragNames, SEXP stutterIndex, SEXP repAdjust, SEXP detectionThresh)
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
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter to double
	SEXP MEANS = PROTECT(duplicate(meanS));
	double const * const meanS_ptr     = REAL(MEANS);
	double means = meanS_ptr[0];

	// convert stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];


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
	SEXP stutterindex = PROTECT(duplicate(stutterIndex));
	double * stutterindex_ptr     = REAL(stutterindex);
	std::vector<double> fragVecN, fragVecL,indexStutter;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(myRound(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		indexStutter.push_back(stutterindex_ptr[i]);
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

	double cdfArg = detectThresh/scaleDouble;
	double pdfArg = 1/scaleDouble; 

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	# pragma omp parallel for schedule(dynamic)
	for(int x=0; x<nCombs; ++x)
		{
		outDouble[x] = singleGenotypeSD(genotypeArrayVec, csp, DNAcontVec, gradients, 
		meand,means, 
		repadjust, scaleDouble, degVec, fragVecL, fragVecN, indexStutter, x, nGen, nCSP, nCont, nFrag, cdfArg, pdfArg);
		//outDouble = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterMeanDouble, stutterAdjustDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, detectThresh, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif

	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    UNPROTECT(14);
	return result;

	}

// get probability of all genotype combinations, taking into account Single and Over stutter
SEXP probabilityPeaksSO(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                    SEXP DNAcont, SEXP gradientS, SEXP meanO,
                    SEXP meanS, 
                    SEXP scale, SEXP degradation, SEXP fragLengths, 
                    SEXP fragNames, SEXP stutterIndex, SEXP repAdjust, SEXP detectionThresh)
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
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter to double
	SEXP MEANS = PROTECT(duplicate(meanS));
	double const * const meanS_ptr     = REAL(MEANS);
	double means = meanS_ptr[0];

	// convert overstutterRate to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];


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
	SEXP stutterindex = PROTECT(duplicate(stutterIndex));
	double * stutterindex_ptr     = REAL(stutterindex);
	std::vector<double> fragVecN, fragVecL,indexStutter;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(myRound(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		indexStutter.push_back(stutterindex_ptr[i]);
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

	double cdfArg = detectThresh/scaleDouble;
	double pdfArg = 1/scaleDouble; 

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	# pragma omp parallel for schedule(dynamic)
	for(int x=0; x<nCombs; ++x)
		{
		outDouble[x] = singleGenotypeSO(genotypeArrayVec, csp, DNAcontVec, gradients, 
		meano,means, 
		repadjust, scaleDouble, degVec, fragVecL, fragVecN, indexStutter, x, nGen, nCSP, nCont, nFrag, cdfArg, pdfArg);
		//outDouble = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterMeanDouble, stutterAdjustDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, detectThresh, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif

	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    UNPROTECT(14);
	return result;

	}

// get probability of all genotype combinations, taking into account Single, Double and Over stutter
SEXP probabilityPeaksSDO(SEXP genotypeArray, SEXP alleles, SEXP heights, 
                    SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO,
                    SEXP meanS, 
                    SEXP scale, SEXP degradation, SEXP fragLengths, 
                    SEXP fragNames, SEXP stutterIndex, SEXP repAdjust, SEXP detectionThresh)
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
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter to double
	SEXP MEANS = PROTECT(duplicate(meanS));
	double const * const meanS_ptr     = REAL(MEANS);
	double means = meanS_ptr[0];

	// convert stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];

	// convert stutter to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];


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
	SEXP stutterindex = PROTECT(duplicate(stutterIndex));
	double * stutterindex_ptr     = REAL(stutterindex);
	std::vector<double> fragVecN, fragVecL,indexStutter;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(myRound(fragNvec_ptr[i]*10.0f)/10.0f);
		fragVecL.push_back(fragLvec_ptr[i]);
		indexStutter.push_back(stutterindex_ptr[i]);
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

	double cdfArg = detectThresh/scaleDouble;
	double pdfArg = 1/scaleDouble; 

	//Rprintf("Before Main Loop");
	// Loop over genotype combinations
	# pragma omp parallel for schedule(dynamic)
	for(int x=0; x<nCombs; ++x)
		{
		outDouble[x] = singleGenotypeSDO(genotypeArrayVec, csp, DNAcontVec, gradients, 
		meand,meano,means, 
		repadjust, scaleDouble, degVec, fragVecL, fragVecN, indexStutter, x, nGen, nCSP, nCont, nFrag, cdfArg, pdfArg);
		//outDouble = singleGenotype(genotypeArrayVec, csp, DNAcontVec, stutterMeanDouble, stutterAdjustDouble, scaleDouble, degVec, fragVecL, fragVecN, repadjust, detectThresh, x, nGen, nCSP, nCont, nFrag);
		}
	//Rprintf("After main loop");
	# ifdef OPENMP_STACK
	//    R_CStackLimit = oldstack;
	# endif

	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    UNPROTECT(15);
	return result;

	}


/*
// get probability of all genotype combinations, taking into account Single, Double and Over stutter
SEXP testCDF(SEXP S, SEXP Z)
    {
    // convert S to double
	SEXP sDuplicate = PROTECT(duplicate(S));
	double const * const S_ptr     = REAL(sDuplicate);
	double s = S_ptr[0];
	// convert Z to double
	SEXP zDuplicate = PROTECT(duplicate(Z));
	double const * const Z_ptr     = REAL(zDuplicate);
	double z = Z_ptr[0];
	// get cdf
	result = ln_kf_gammap(s,z)
	// Make and return output object
	SEXP out = PROTECT(allocVector(REALSXP, 1));
  	double       * const out_ptr  = REAL(out);
	out_ptr[0] = result;
    UNPROTECT(3);
	return(out);
    }
*/








