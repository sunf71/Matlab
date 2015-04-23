/*
 * MATLAB Compiler: 4.14 (R2010b)
 * Date: Wed Dec 03 14:18:45 2014
 * Arguments: "-B" "macro_default" "-W" "lib:libAnn" "-T" "link:lib" "CSH_nn.m" 
 */

#include <stdio.h>
#define EXPORTING_libAnn 1
#include "libAnn.h"

static HMCRINSTANCE _mcr_inst = NULL;


#if defined( _MSC_VER) || defined(__BORLANDC__) || defined(__WATCOMC__) || defined(__LCC__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#include <windows.h>

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libAnn_C_API
#define LIB_libAnn_C_API /* No special import/export declaration */
#endif

LIB_libAnn_C_API 
bool MW_CALL_CONV libAnnInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!GetModuleFileName(GetModuleHandle("libAnn"), path_to_dll, _MAX_PATH))
    return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll, 
                                    1728348);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(   &_mcr_inst,
                                                                error_handler, 
                                                                print_handler,
                                                                ctfStream, 
                                                                1728348);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
  return true;
}

LIB_libAnn_C_API 
bool MW_CALL_CONV libAnnInitialize(void)
{
  return libAnnInitializeWithHandlers(mclDefaultErrorHandler, mclDefaultPrintHandler);
}

LIB_libAnn_C_API 
void MW_CALL_CONV libAnnTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

LIB_libAnn_C_API 
long MW_CALL_CONV libAnnGetMcrID() 
{
  return mclGetID(_mcr_inst);
}

LIB_libAnn_C_API 
void MW_CALL_CONV libAnnPrintStackTrace(void) 
{
  char** stackTrace;
  int stackDepth = mclGetStackTrace(&stackTrace);
  int i;
  for(i=0; i<stackDepth; i++)
  {
    mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
    mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
  }
  mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_libAnn_C_API 
bool MW_CALL_CONV mlxCSH_nn(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
  return mclFeval(_mcr_inst, "CSH_nn", nlhs, plhs, nrhs, prhs);
}

LIB_libAnn_C_API 
bool MW_CALL_CONV mlfCSH_nn(int nargout, mxArray** CSH_ann, mxArray** CSH_bnn, mxArray* 
                            A, mxArray* B, mxArray* width, mxArray* iterations, mxArray* 
                            k, mxArray* calcBnn, mxArray* bMask, mxArray* patch_mode, 
                            mxArray* patch_params)
{
  return mclMlfFeval(_mcr_inst, "CSH_nn", nargout, 2, 9, CSH_ann, CSH_bnn, A, B, width, iterations, k, calcBnn, bMask, patch_mode, patch_params);
}
