#ifndef __LIKELTD_GAMMADIST_H_INCLUDED__  
#define __LIKELTD_GAMMADIST_H_INCLUDED__  
 

#include "config.h"

#include <R.h>
#include <Rdefines.h>
#include <vector>
#include <string>
 
#ifdef __cplusplus
extern "C" {
#endif
    //! \brief Error handling function
    //void nrerror(char error_text);

    //! \brief Returns the log gamma function
    double new_gammln(double xx);

	//! \brief Returns the incomplete gamma function P(a, x).
    double gammp(double a, double x);

    //! \brief Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
    double gammq(double a, double x);

    //! \brief Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
    void gser(double *gamser, double a, double x, double *gln);

    //! \brief Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
    void gcf(double *gammcf, double a, double x, double *gln);


#ifdef __cplusplus
}
#endif
#endif
