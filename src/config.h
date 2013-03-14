#ifndef LIKELTD_CONFIG_H
#define LIKELTD_CONFIG_H

#if defined(__APPLE__) || defined(__MACH__)
#  ifdef _OPENMP
#    define OPENMP_STACK
#  endif 
#elif defined(__linux__) || defined(__unix) || defined(__linux)
#  ifdef _OPENMP
#    define OPENMP_STACK
#  endif
#endif

#endif
