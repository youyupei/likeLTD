// Functions for gamma distribution from Numerical Recipes
#include "config.h"
#include "openmp.h"
#include "maximizePeaks.h"
#include "gammaDist.h"


#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <R_ext/Error.h>

#include <math.h>
#define ITMAX 1.0e6  // maximum allowed number of iterations
#define EPS 3.0e-7  // relative accuracy
#define FPMIN 1.0e-30  // Number near smallest representable floating point number 

// Numerical Recipes standard error handler
void nrerror(char error_text[])
    {
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
    }

// Returns the value ln[Γ(xx)] for xx > 0.
// Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
// accuracy is good enough.
float new_gammln(float xx)
    {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
    }


//Returns the incomplete gamma function P(a, x).
float gammp(float a, float x)
    {
    void gcf(float *gammcf, float a, float x, float *gln);
    void gser(float *gamser, float a, float x, float *gln);
    void nrerror(char error_text[]);
    float gamser,gammcf,gln;
    if((a>=0.0&&a<1.0e-30)||std::isinf(a)) return kf_gammap(a,x);
    if (x < 0.0 || a < 0.0) 
        {
        fprintf(stderr,"a=%f\nx=%f\n",a,x);
	char message[] = "Invalid arguments in routine gammp";
        nrerror(message);
        }
    if (x < (a+1.0)) { //Use the series representation.
        gser(&gamser,a,x,&gln);
        return gamser;
        } else { //Use the continued fraction representation
        gcf(&gammcf,a,x,&gln);
        return 1.0-gammcf; //and take its complement.
        }
    }

//Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
float gammq(float a, float x)
    {
    void gcf(float *gammcf, float a, float x, float *gln);
    void gser(float *gamser, float a, float x, float *gln);
    void nrerror(char error_text[]);
    float gamser,gammcf,gln;
    if (x < 0.0 || a <= 0.0) 
	{
	char message[] = "Invalid arguments in routine gammq";
	nrerror(message);
	}
    if (x < (a+1.0)) { //Use the series representation
        gser(&gamser,a,x,&gln);
        return 1.0-gamser; //and take its complement.
        } else { //Use the continued fraction representation.
        gcf(&gammcf,a,x,&gln);
        return gammcf;
        }
    }



// Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
// Also returns ln Γ(a) as gln.
void gser(float *gamser, float a, float x, float *gln)
    {
    float new_gammln(float xx);
    void nrerror(char error_text[]);
    int n;
    float sum,del,ap;
    *gln=new_gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) 
		{
		char message[] = "x less than 0 in routine gser";
		nrerror(message);
		}
        *gamser=0.0;
        return;
        } else {
        ap=a;
        del=sum=1.0/a;
        for (n=1;n<=ITMAX;n++) {
            ++ap;
            del *= x/ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS) {
                *gamser=sum*exp(-x+a*log(x)-(*gln));
                return;
                }
            }
        fprintf(stderr,"a=%.30f\nx=%f\n",a,x);
	char message[] = "a too large, ITMAX too small in routine gser";
	nrerror(message);
        return;
        }
    }

//Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
//as gammcf. Also returns lnΓ(a) as gln.
void gcf(float *gammcf, float a, float x, float *gln)
    {
    float new_gammln(float xx);
    void nrerror(char error_text[]);
    int i;
    float an,b,c,d,del,h;
    *gln=new_gammln(a);
    b=x+1.0-a; //Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 = 0.
    c=1.0/FPMIN;
    d=1.0/b;
    h=d;
    for (i=1;i<=ITMAX;i++) { //Iterate to convergence.
        an = -i*(i-a);
        b += 2.0;
        d=an*d+b;
        if (fabs(d) < FPMIN) d=FPMIN;
        c=b+an/c;
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
        }
    if (i > ITMAX) 
	{
	char message[] = "a too large, ITMAX too small in gcf";
	nrerror(message);
	}
    *gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
    }



















