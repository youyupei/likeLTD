
#include "config.h"
#include "openmp.h"
#include "maximizePeaks.h"
#include "gammaDist.h"

#include <cmath>
#include <string>
#include <vector>
#include <iterator>
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

//inline double roundOneDP(double d) 
//    {
//    return myRound(d*10.0f)/10.0f;
//    }

inline double roundOneDP(double d) 
    {
    return floor(d*10+0.5)/10;
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



// regularized lower incomplete gamma function, by series expansion
// see: http://barricklab.org/hg/breseq/file/9ab26de12a0b/extern/samtools-0.1.15/bcftools/kfunc.c
#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290
// regularized lower incomplete gamma function, by series expansion
double _kf_gammap(double s, double z)
{
	double sum, x;
 	int k;
 	for (k = 1, sum = x = 1.; k < 100; ++k) {
 		sum += (x *= z / (s + k));
 		if (x / sum < KF_GAMMA_EPS) break;
 	}
 	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
 }
 // regularized upper incomplete gamma function, by continued fraction
double _kf_gammaq(double s, double z)
 {
 	int j;
 	double C, D, f;
 	f = 1. + z - s; C = f; D = 0.;
 	// Modified Lentz's algorithm for computing continued fraction
 	// See Numerical Recipes in C, 2nd edition, section 5.2
 	for (j = 1; j < 100; ++j) {
 		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
 		D = b + a * D;
 		if (D < KF_TINY) D = KF_TINY;
 		C = b + a / C;
 		if (C < KF_TINY) C = KF_TINY;
 		D = 1. / D;
 		d = C * D;
 		f *= d;
 		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
 	}
 	return exp(s * log(z) - z - kf_lgamma(s) - log(f));
 }

double kf_gammap(double s, double z)
 {
 	return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
 }

double ln_kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return s * log(z) - z - kf_lgamma(s + 1.) + log(sum);
}


SEXP testCDF(SEXP S, SEXP Z)
    {
    // convert S to double
	double const * const S_ptr     = REAL(S);
	double s = S_ptr[0];
	// convert Z to double
	double const * const Z_ptr     = REAL(Z);
	double z = Z_ptr[0];
	// get cdf
	double result = kf_gammap(s,z);
	// Make and return output object
	SEXP out = PROTECT(allocVector(REALSXP, 1));
  	double       * const out_ptr  = REAL(out);
	out_ptr[0] = result;
    UNPROTECT(1);
	return(out);
    }

SEXP testNewCDF(SEXP S, SEXP Z)
    {
    // convert S to double
	double const * const S_ptr     = REAL(S);
	double s = S_ptr[0];
	// convert Z to double
	double const * const Z_ptr     = REAL(Z);
	double z = Z_ptr[0];
	// get cdf
	double result = gammp(s, z);
	// Make and return output object
	SEXP out = PROTECT(allocVector(REALSXP, 1));
  	double       * const out_ptr  = REAL(out);
	out_ptr[0] = result;
    UNPROTECT(1);
	return(out);
    }

SEXP testPDF(SEXP X, SEXP A, SEXP B)
    {
    // convert X to double
	SEXP xDuplicate = PROTECT(duplicate(X));
	double const * const X_ptr     = REAL(xDuplicate);
	double x = X_ptr[0];
	// convert A to double
	SEXP aDuplicate = PROTECT(duplicate(A));
	double const * const A_ptr     = REAL(aDuplicate);
	double a = A_ptr[0];
	// convert B to double
	SEXP bDuplicate = PROTECT(duplicate(B));
	double const * const B_ptr     = REAL(bDuplicate);
	double b = B_ptr[0];
	// get cdf
	double result = gammalog(x,a,b);
	// Make and return output object
	SEXP out = PROTECT(allocVector(REALSXP, 1));
  	double       * const out_ptr  = REAL(out);
	out_ptr[0] = result;
    UNPROTECT(4);
	return(out);
    }












SEXP getProbabilitiesSDO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals)
	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert double stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];

	// convert over stutter to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];

	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

	// convert fragLengths, fragNames, LUSvals to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
	// get probability for each genotype combination
	# pragma omp parallel for //schedule(dynamic)
	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*4,0);

        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// double stutter positions
            		doubleStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-2.0);
    	    		// over stutter positions
            		overStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]+1.0);
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
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));

        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double doubleStutterDose, overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);
				
					if(diff<0.0001) 
				    		{
				   	 	matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				//Rprintf("%f\n",degVec[i]);
				// compute linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
            			// stutter doses
				stutterDose = tmpDose * stutterRate;
				doubleStutterDose = tmpDose * meand;
				overStutterDose = tmpDose * meano;
				nonstutterDose = tmpDose * (1-(stutterRate+meand+meano));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// double stutter dose
			tmpMu.genotype = doubleStutterVec[x];
			tmpMu.dose = doubleStutterDose;
			muSd[x] = tmpMu;
			// over stutter dose
			tmpMu.genotype = overStutterVec[x];
			tmpMu.dose = overStutterDose;
			muSo[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(std::fabs(muA[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muA[j].dose;
				if(std::fabs(muS[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muS[j].dose;
				if(std::fabs(muSd[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muSd[j].dose;
				if(std::fabs(muSo[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muSo[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
					gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(17);
	return result;
    	}


SEXP getProbabilitiesSO(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals)
	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));
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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert over stutter to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];

	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

	// convert fragLengths, fragNames, LUSvals to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
	// get probability for each genotype combination
	# pragma omp parallel for //schedule(dynamic)
	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*3,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// over stutter positions
            		overStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]+1.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			// over stutter
    			allPosVec[y+(nGen*2)] = overStutterVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				// linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
				// stutter doses
				stutterDose = tmpDose * stutterRate;
				overStutterDose = tmpDose * meano;
				nonstutterDose = tmpDose * (1-(stutterRate+meano));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// over stutter dose
			tmpMu.genotype = overStutterVec[x];
			tmpMu.dose = overStutterDose;
			muSo[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				if(muSo[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muSo[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(16);
	return result;
    	}


SEXP getProbabilitiesSD(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals)
	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert double stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];

	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), allPosVec(nGen*3,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// double stutter positions
            		doubleStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-2.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			// double stutter
    			allPosVec[y+(nGen*2)] = doubleStutterVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double doubleStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);	
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
            			// get linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
            			// stutter doses
				stutterDose = tmpDose * stutterRate;
				doubleStutterDose = tmpDose * meand;
				nonstutterDose = tmpDose * (1-(stutterRate+meand));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// double stutter dose
			tmpMu.genotype = doubleStutterVec[x];
			tmpMu.dose = doubleStutterDose;
			muSd[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				if(muSd[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muSd[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}     
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(16);
	return result;
    	}


SEXP getProbabilitiesS(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals)
    	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
    	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), allPosVec(nGen*2,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);	
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				// linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
				// stutter doses
				stutterDose = tmpDose * stutterRate;
				nonstutterDose = tmpDose * (1-stutterRate);
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
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


SEXP getProbabilitiesSDO_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg)
	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert double stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];

	// convert over stutter to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];

	// convert dropin to double
	SEXP DROPIN = PROTECT(duplicate(dropin));
	double const * const dropin_ptr     = REAL(DROPIN);
	double Dropin = dropin_ptr[0];

	// convert dropinDeg to double
	SEXP DROPINDEG = PROTECT(duplicate(dropinDeg));
	double const * const dropindeg_ptr     = REAL(DROPINDEG);
	double Dropindeg = dropindeg_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
    	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals and fragProbs to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP fragPvec = PROTECT(duplicate(fragProbs));
	double * fragPvec_ptr     = REAL(fragPvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		fragVecP.push_back(fragPvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*4,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// double stutter positions
            		doubleStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-2.0);
    	    		// over stutter positions
            		overStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]+1.0);
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
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double doubleStutterDose, overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				//Rprintf("%f\n",degVec[i]);
				// compute linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
            			// stutter doses
				stutterDose = tmpDose * stutterRate;
				doubleStutterDose = tmpDose * meand;
				overStutterDose = tmpDose * meano;
				nonstutterDose = tmpDose * (1-(stutterRate+meand+meano));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// double stutter dose
			tmpMu.genotype = doubleStutterVec[x];
			tmpMu.dose = doubleStutterDose;
			muSd[x] = tmpMu;
			// over stutter dose
			tmpMu.genotype = overStutterVec[x];
			tmpMu.dose = overStutterDose;
			muSo[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(std::fabs(muA[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muA[j].dose;
				if(std::fabs(muS[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muS[j].dose;
				if(std::fabs(muSd[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muSd[j].dose;
				if(std::fabs(muSo[j].genotype-allPosVec[x])<0.0001) tmpDose2 = tmpDose2+muSo[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
		// add dropin doses to dose array
		for(int j=0; j<fragVecN.size(); j++)
			{
			if(!(fragVecN[j]<-1&&fragVecN[j]>-100))
		    		{
		    		int matchIndex=0;
		    		for(int k=0; k<dbVals.size(); k++)
			    		{
			    		double diff = std::fabs(dbVals[k]-fragVecN[j]); 
                			if(diff<0.0001)
	                			{
        	        			matchIndex = k;
        	        			break;
        	        			}
			    		}
		    		doseArray[matchIndex][i] += fragVecP[j]*Dropin*pow(Dropindeg,-fragVecL[j]);
		    		}
			}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif		
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(20);
	return result;
    	}


SEXP getProbabilitiesSO_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanO, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg)
    	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert over stutter to double
	SEXP MEANO = PROTECT(duplicate(meanO));
	double const * const meanO_ptr     = REAL(MEANO);
	double meano = meanO_ptr[0];

	// convert dropin to double
	SEXP DROPIN = PROTECT(duplicate(dropin));
	double const * const dropin_ptr     = REAL(DROPIN);
	double Dropin = dropin_ptr[0];

	// convert dropinDeg to double
	SEXP DROPINDEG = PROTECT(duplicate(dropinDeg));
	double const * const dropindeg_ptr     = REAL(DROPINDEG);
	double Dropindeg = dropindeg_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
    	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals and fragProbs to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP fragPvec = PROTECT(duplicate(fragProbs));
	double * fragPvec_ptr     = REAL(fragPvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		fragVecP.push_back(fragPvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), overStutterVec(nGen,0), allPosVec(nGen*3,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// over stutter positions
            		overStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]+1.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			// over stutter
    			allPosVec[y+(nGen*2)] = overStutterVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSo(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				// linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
				// stutter doses
				stutterDose = tmpDose * stutterRate;
				overStutterDose = tmpDose * meano;
				nonstutterDose = tmpDose * (1-(stutterRate+meano));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// over stutter dose
			tmpMu.genotype = overStutterVec[x];
			tmpMu.dose = overStutterDose;
			muSo[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				if(muSo[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muSo[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
		// add dropin doses to dose array
		for(int j=0; j<fragVecN.size(); j++)
			{
			if(!(fragVecN[j]<-1&&fragVecN[j]>-100))
		    		{
		    		int matchIndex=0;
		    		for(int k=0; k<dbVals.size(); k++)
			    		{
			    		double diff = std::fabs(dbVals[k]-fragVecN[j]); 
                			if(diff<0.0001)
	                			{
        	        			matchIndex = k;
        	        			break;
        	        			}
			    		}
		    		//Rprintf("%d\t",matchIndex);
		    		doseArray[matchIndex][i] += fragVecP[j]*Dropin*pow(Dropindeg,-fragVecL[j]);
		    		}
			}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(19);
	return result;
    	}

SEXP getProbabilitiesSD_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP meanD, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg)
    	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert double stutter to double
	SEXP MEAND = PROTECT(duplicate(meanD));
	double const * const meanD_ptr     = REAL(MEAND);
	double meand = meanD_ptr[0];

	// convert dropin to double
	SEXP DROPIN = PROTECT(duplicate(dropin));
	double const * const dropin_ptr     = REAL(DROPIN);
	double Dropin = dropin_ptr[0];

	// convert dropinDeg to double
	SEXP DROPINDEG = PROTECT(duplicate(dropinDeg));
	double const * const dropindeg_ptr     = REAL(DROPINDEG);
	double Dropindeg = dropindeg_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
    	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals and fragProbs to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP fragPvec = PROTECT(duplicate(fragProbs));
	double * fragPvec_ptr     = REAL(fragPvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		fragVecP.push_back(fragPvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), doubleStutterVec(nGen,0), allPosVec(nGen*3,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    	    		// double stutter positions
            		doubleStutterVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-2.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			// double stutter
    			allPosVec[y+(nGen*2)] = doubleStutterVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu),muSd(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double doubleStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
            			// get linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
            			// stutter doses
				stutterDose = tmpDose * stutterRate;
				doubleStutterDose = tmpDose * meand;
				nonstutterDose = tmpDose * (1-(stutterRate+meand));
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			// double stutter dose
			tmpMu.genotype = doubleStutterVec[x];
			tmpMu.dose = doubleStutterDose;
			muSd[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				if(muSd[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muSd[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{		       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
		// add dropin doses to dose array
		for(int j=0; j<fragVecN.size(); j++)
			{
			if(!(fragVecN[j]<-1&&fragVecN[j]>-100))
		    		{
		    		int matchIndex=0;
		    		for(int k=0; k<dbVals.size(); k++)
			    		{
			    		double diff = std::fabs(dbVals[k]-fragVecN[j]); 
                			if(diff<0.0001)
	                			{
        	        			matchIndex = k;
        	        			break;
        	        			}
			    		}
		    		//Rprintf("%d\t",matchIndex);
		   		doseArray[matchIndex][i] += fragVecP[j]*Dropin*pow(Dropindeg,-fragVecL[j]);
		    		}
			}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}      
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(19);
	return result;
    	}


SEXP getProbabilitiesS_dropin(SEXP genotypeArray, SEXP DNAcont, SEXP gradientS, SEXP interceptS, SEXP degradation, SEXP fragLengths, SEXP fragNames, SEXP LUSvals, SEXP alleles, SEXP heights, SEXP repAdjust, SEXP scale, SEXP detectionThresh, SEXP databaseVals,SEXP fragProbs,SEXP dropin,SEXP dropinDeg)
    	{
    	# ifdef OPENMP_STACK
	//    uintptr_t const oldstack = R_CStackLimit;
	//    R_CStackLimit = (uintptr_t) - 1;
	# endif
	// get various params
	int const nCombs = INTEGER(GET_DIM(genotypeArray))[1];
	int const nGen = INTEGER(GET_DIM(genotypeArray))[0];
	int const nCont = length(DNAcont);
	int const nFrag = length(fragLengths);
	int const nDeg = length(degradation);
        int const nDat = length(databaseVals);
	// create some objects
    	std::vector<std::vector<double> > doseArray(nDat, std::vector<double>(nCombs));

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

	// convert stutter gradient to double
	SEXP GRADIENTS = PROTECT(duplicate(gradientS));
	double const * const gradientS_ptr     = REAL(GRADIENTS);
	double gradients = gradientS_ptr[0];

	// convert stutter intercept to double
	SEXP INTERCEPTS = PROTECT(duplicate(interceptS));
	double const * const interceptS_ptr     = REAL(INTERCEPTS);
	double intercepts = interceptS_ptr[0];

	// convert dropin to double
	SEXP DROPIN = PROTECT(duplicate(dropin));
	double const * const dropin_ptr     = REAL(DROPIN);
	double Dropin = dropin_ptr[0];

	// convert dropinDeg to double
	SEXP DROPINDEG = PROTECT(duplicate(dropinDeg));
	double const * const dropindeg_ptr     = REAL(DROPINDEG);
	double Dropindeg = dropindeg_ptr[0];

    	// convert degradation to vector
	SEXP DEGRADATION = PROTECT(duplicate(degradation));
	double const * const deg_ptr     = REAL(DEGRADATION);
	std::vector<double> degVec;
	for(int i=0; i<nDeg; ++i)
		{
		degVec.push_back(deg_ptr[i]);
		}
		
    	// convert databaseVals to vector
	SEXP DATABASEVALS = PROTECT(duplicate(databaseVals));
	double const * const dbvals_ptr     = REAL(DATABASEVALS);
	std::vector<double> dbVals;
	for(int i=0; i<nDat; ++i)
		{
		dbVals.push_back(dbvals_ptr[i]);
		}

	// get more parameters
	int const nRep = length(repAdjust);
	// create more objects
	std::vector<double> tmpVec;
	std::vector<double> outDouble;
	outDouble.assign(nCombs,1);

	// convert alleles to vector of vectors (equivalent to list)
	SEXP ALLELES = PROTECT(duplicate(alleles));
	std::vector<std::vector<double> >  allelesVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(ALLELES,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(ALLELES,r))[i]);
		    	}	
		allelesVec.push_back(tmpVec);    
		}

	// convert heights to vector of vectors (equivalent to list)
	SEXP HEIGHTS = PROTECT(duplicate(heights));
	std::vector<std::vector<double> >  heightsVec;
	for(int r=0; r<nRep; r++)
	    	{
		tmpVec.clear();
	    	for(int i=0; i<length(VECTOR_ELT(HEIGHTS,r)); ++i)
		    	{
		    	tmpVec.push_back(REAL(VECTOR_ELT(HEIGHTS,r))[i]);
		    	}	
		heightsVec.push_back(tmpVec);    
		}

	// convert repAdjust to vector
	SEXP REPADJUST = PROTECT(duplicate(repAdjust));
	double const * const repadjust_ptr     = REAL(REPADJUST);
	std::vector<double> repadjustVec(nRep,0);
	for(int i=0; i<nRep; ++i)
		{
		repadjustVec[i] = repadjust_ptr[i];
		}	

	// convert scale to double
	SEXP SCALE = PROTECT(duplicate(scale));
	double const * const scale_ptr     = REAL(SCALE);
	double scaleDouble = scale_ptr[0];

    	// convert detectionThresh to double
	SEXP DETECTIONTHRESH = PROTECT(duplicate(detectionThresh));
	double const * const detect_ptr     = REAL(DETECTIONTHRESH);
	double detectDouble = detect_ptr[0];

    	// convert fragLengths, fragNames, LUSvals and fragProbs to vector
	SEXP fragNvec = PROTECT(duplicate(fragNames));
	double * fragNvec_ptr     = REAL(fragNvec);
	SEXP fragLvec = PROTECT(duplicate(fragLengths));
	double * fragLvec_ptr     = REAL(fragLvec);
	SEXP fragPvec = PROTECT(duplicate(fragProbs));
	double * fragPvec_ptr     = REAL(fragPvec);
	SEXP lusvals = PROTECT(duplicate(LUSvals));
	double * lusvals_ptr     = REAL(lusvals);
	std::vector<double> fragVecN, fragVecL, fragVecP,lusVals;
	for(int i=0; i<nFrag; ++i)
		{
		fragVecN.push_back(fragNvec_ptr[i]);
		fragVecL.push_back(fragLvec_ptr[i]);
		fragVecP.push_back(fragPvec_ptr[i]);
		lusVals.push_back(lusvals_ptr[i]);	
		}
	// get argument for cdf
	double cdfArg = detectDouble/scaleDouble;
    	// get probability for each genotype combination
    	# pragma omp parallel for //schedule(dynamic)
    	for(int i=0; i<nCombs; i++)
        	{
		// loop over genotype combinations
        	std::vector<double> genotypeVec(nGen,0), stutterPosVec(nGen,0), allPosVec(nGen*2,0);
        	// Loop over members of genotype
	    	for(int y=0; y<nGen; ++y)
    			{
    	    		// round genotypes
    			genotypeVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]);
    			// get stutter positions
		    	stutterPosVec[y] = roundOneDP(genotypeArrayVec[(i*nGen)+y]-1.0);
    			// get all positions, while rounding to 1dp
    			allPosVec[y] = genotypeVec[y];
    			allPosVec[y+nGen] = stutterPosVec[y];
    			}
    		// sort and unique allPos
	    	std::sort(allPosVec.begin(),allPosVec.end());
    		std::vector<double>::iterator itFlt = std::unique(allPosVec.begin(),allPosVec.end());
    		allPosVec.resize(std::distance(allPosVec.begin(),itFlt));
        	// get doses ignoring repAdjust
		double tmpDose, fragSub, stuttIndSub, stutterDose, nonstutterDose;
		genoStruct tmpMu;
		std::vector<genoStruct> muA(genotypeVec.size(),tmpMu),muS(genotypeVec.size(),tmpMu);
		int matchIndex;
		double diff;
		double overStutterDose, stutterRate;
		// effective dose for stutter and allelic
		for(int x=0; x<nGen; ++x)
			{
			// Do not compute dose again if homozygote
			if((x>0)&&(x%2!=0)&&(genotypeVec[x]==genotypeVec[x-1]))
				{
			
				} else {
				// get fragLengths
				matchIndex = 0;
				for(unsigned int q=0; q<nFrag; ++q)
					{
					diff = std::fabs(genotypeVec[x]-fragVecN[q]);	
					if(diff<0.0001) 
				    		{
				    		matchIndex = q;
	                			break;
				    		}
					}
				fragSub = fragVecL[matchIndex];
				stuttIndSub = lusVals[matchIndex];
				// compute effective dose
				tmpDose = DNAcontVec[x]*std::pow(degVec[x],-fragSub);
				// linear stutter rate
				stutterRate = intercepts+(gradients*stuttIndSub);
				// stutter doses
				stutterDose = tmpDose * stutterRate;
				nonstutterDose = tmpDose * (1-stutterRate);
				}
			// stutter adjusted effective dose
			// non-stutter dose
			tmpMu.genotype = genotypeVec[x];
			tmpMu.dose = nonstutterDose;
			muA[x] = tmpMu;
			// stutter dose
			tmpMu.genotype = stutterPosVec[x];
			tmpMu.dose = stutterDose;
			muS[x] = tmpMu;
			}
		genoStruct tmpMu2; 
		std::vector<genoStruct> gammaMuVec(allPosVec.size(),tmpMu2);
		double tmpDose2;
		// combine doses at each allelic position
		for(unsigned int x=0; x<allPosVec.size(); ++x)
			{
			tmpDose2 = 0;
			for(unsigned int j=0; j<muA.size(); ++j)
				{
				if(muA[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muA[j].dose;
				if(muS[j].genotype==allPosVec[x]) tmpDose2 = tmpDose2+muS[j].dose;
				}
			tmpMu2.genotype = allPosVec[x];
			tmpMu2.dose = tmpDose2;
			gammaMuVec[x] = tmpMu2;
			}
        	// slot doses into dose array
        	for(int j=0; j<gammaMuVec.size(); j++)
            		{       
            		int matchIndex = 0;         
            		for(int k=0; k<dbVals.size(); k++)
                		{
                		double diff = std::fabs(dbVals[k]-gammaMuVec[j].genotype); 
                		if(diff<0.0001)
                    			{
                    			matchIndex = k;
                    			break;
                    			}
                		}
            		doseArray[matchIndex][i] += gammaMuVec[j].dose;
            		}
		// add dropin doses to dose array
		for(int j=0; j<fragVecN.size(); j++)
			{
			if(!(fragVecN[j]<-1&&fragVecN[j]>-100))
		    		{
		    		int matchIndex=0;
		    		for(int k=0; k<dbVals.size(); k++)
			    		{
			    		double diff = std::fabs(dbVals[k]-fragVecN[j]); 
                			if(diff<0.0001)
	                			{
        	        			matchIndex = k;
        	        			break;
        	        			}
			    		}
		    		//Rprintf("%d\t",matchIndex);
		    		doseArray[matchIndex][i] += fragVecP[j]*Dropin*pow(Dropindeg,-fragVecL[j]);
		    		}
			}
    		// get probabilities
		// loop over database alleles
		for(int j=0; j<nDat; j++)
			{
			// loop over replicates
			for(int k=0; k<nRep; k++)
		    		{
	        		int matchIndex = 0;   
	        		bool matchFlag = false; 
		    		// only do something if some hypothesised dose
		    		if(doseArray[j][i]!='\0'&&doseArray[j][i]!=0)
		        		{
		        		// only do anything if doseArray is not 0
		        		// which allele are we looking at?
		        		int matchIndex = 0;   
		        		bool matchFlag = false;     
        				for(int l=0; l<allelesVec[k].size(); l++)
        	                		{
        	                		double diff = std::fabs(dbVals[j]-allelesVec[k][l]);  
        	                		if(diff<0.0001)
        	                    			{           
        	                    			matchIndex = l;
        	                    			matchFlag = true;
        	                    			break;
        	                    			}
        	                		}
					if(matchFlag==false)
		    				{
						double prob = gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,cdfArg);
							}
                    				// dropout dose
                    				outDouble[i] = outDouble[i] * prob;
		    				} else {
						double prob = (gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
						gammp((doseArray[j][i]*repadjustVec[k])/scaleDouble,
						(heightsVec[k][matchIndex]-0.5)/scaleDouble));
						if(std::fabs(prob-999)<0.0001) // if error in gammp
							{
							prob = (kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]+0.5)/scaleDouble)-
							kf_gammap((doseArray[j][i]*repadjustVec[k])/scaleDouble,
							(heightsVec[k][matchIndex]-0.5)/scaleDouble));
							}
                    				// non-dropout dose
                    				outDouble[i] = outDouble[i] * prob;
			            		}   
			        	}
			    	}
			}
		//Rprintf("After main loop");
	    	# ifdef OPENMP_STACK
	    	//    R_CStackLimit = oldstack;
	    	# endif	
		}
	// Make and return output object
	SEXP result = PROTECT(allocVector(REALSXP, nCombs));
  	double       * const out_ptr  = REAL(result);
	for(int i=0; i<nCombs; ++i)
		{
		out_ptr[i] = outDouble[i];
		}
    	UNPROTECT(18);
	return result;
    	}




