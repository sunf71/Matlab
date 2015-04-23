/*
 * MATLAB Compiler: 4.14 (R2010b)
 * Date: Wed Dec 03 14:18:45 2014
 * Arguments: "-B" "macro_default" "-W" "lib:libAnn" "-T" "link:lib" "CSH_nn.m" 
 */

#ifndef __libAnn_h
#define __libAnn_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libAnn
#define PUBLIC_libAnn_C_API __global
#else
#define PUBLIC_libAnn_C_API /* No import statement needed. */
#endif

#define LIB_libAnn_C_API PUBLIC_libAnn_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libAnn
#define PUBLIC_libAnn_C_API __declspec(dllexport)
#else
#define PUBLIC_libAnn_C_API __declspec(dllimport)
#endif

#define LIB_libAnn_C_API PUBLIC_libAnn_C_API


#else

#define LIB_libAnn_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libAnn_C_API 
#define LIB_libAnn_C_API /* No special import/export declaration */
#endif

extern LIB_libAnn_C_API 
bool MW_CALL_CONV libAnnInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libAnn_C_API 
bool MW_CALL_CONV libAnnInitialize(void);

extern LIB_libAnn_C_API 
void MW_CALL_CONV libAnnTerminate(void);



extern LIB_libAnn_C_API 
void MW_CALL_CONV libAnnPrintStackTrace(void);

extern LIB_libAnn_C_API 
bool MW_CALL_CONV mlxCSH_nn(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

extern LIB_libAnn_C_API 
long MW_CALL_CONV libAnnGetMcrID();



extern LIB_libAnn_C_API bool MW_CALL_CONV mlfCSH_nn(int nargout, mxArray** CSH_ann, mxArray** CSH_bnn, mxArray* A, mxArray* B, mxArray* width, mxArray* iterations, mxArray* k, mxArray* calcBnn, mxArray* bMask, mxArray* patch_mode, mxArray* patch_params);

#ifdef __cplusplus
}
#endif
#endif
