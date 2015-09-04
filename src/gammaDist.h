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
    void nrerror(char error_text);

    //! \brief Returns the log gamma function
    float new_gammln(float xx);

	//! \brief Returns the incomplete gamma function P(a, x).
    float gammp(float a, float x);

    //! \brief Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
    float gammq(float a, float x);

    //! \brief Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
    void gser(float *gamser, float a, float x, float *gln);

    //! \brief Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
    void gcf(float *gammcf, float a, float x, float *gln);


#ifdef __cplusplus
}
#endif
#endif