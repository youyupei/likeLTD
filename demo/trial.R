code = "SEXP res;
        int const nrow = INTEGER(GET_DIM(x))[0];
        int const ncol = INTEGER(GET_DIM(x))[1];
        PROTECT(res = allocVector(INTSXP, nrow));
        int *inptr = INTEGER(x);
        int * const outptr_first = INTEGER(res);
        int *outptr = outptr_first;
        // Copy first column into output vector.
        for(int i=0; i < nrow; ++i, ++inptr, ++outptr)
          *outptr = (*inptr == NA_INTEGER) ? 1: *inptr;
        // Perform multiplication with other columns
        for(int j=1; j < ncol; ++j) {
          outptr = outptr_first;
          for(int i=0; i < nrow; ++i, ++inptr, ++outptr)
            if(*inptr != NA_INTEGER) *outptr *= *inptr;
        }
        UNPROTECT(1);
        return res;"
prod.matrix.C = cfunction(signature(x="array"), code)
 
prod.matrix.R <- function(x) {
  # Fast row-wise produce.
  #
  # Sadly, this is faster than apply(x, 1, prod)  
  y=x[,1]
  for(i in 2:dim(x)[2])
  y=y*x[,i]
  return(y)
}

library(microbenchmark)
x <- matrix(1:4, 6, 4)
microbenchmark( prod.matrix.C(x),
                prod.matrix.R(x) )

source('R/genetics.R')
source('R/objectives.R')
na = 6
nconfs = 100000
freq = runif(nconfs, 0.1, 0.9)
freqMat = matrix(freq, nrow=nconfs, ncol=na, byrow=F)
unc = rep(FALSE, na)
csp = c(TRUE, rep(FALSE, na-2), TRUE)
res = array(1, nconfs)
zeroAll = matrix(runif(nconfs*na, 0, 1), ncol=nconfs)
dropout = 0.5
