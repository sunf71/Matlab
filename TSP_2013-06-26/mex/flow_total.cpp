// =============================================================================
// == flow_total.cpp
// == --------------------------------------------------------------------------
// == A MEX interface to convert flows for each SP to two images. See m files
// == for usage.
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================

#include "mex.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "debugMEX.h"
#include "helperMEX.h"
#include "matrix.h"

#define NUMARGS 6
#define NUMOUT 2


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, SCALAR ,  "double"); //K
   checkInput(prhs, 1, SCALAR ,  "double"); //xdim
   checkInput(prhs, 2, SCALAR ,  "double"); //ydim
   checkInput(prhs, 3, STRUCT            ); //indices
   checkInput(prhs, 4, VECTOR ,  "double"); //u
   checkInput(prhs, 5, VECTOR ,  "double"); //v


   int K = getInput<double>(prhs[0]);
   int xdim = getInput<double>(prhs[1]);
   int ydim = getInput<double>(prhs[2]);

   arr( arr(double) ) indices_all = allocate_memory<arr(double)>(K);
   arr(double) indices_N = allocate_memory<double>(K);
   for (int k=0; k<K; k++)
   {
      indices_all[k] = getArrayInput<double>(mxGetField(prhs[3],k,"all"));
      indices_N[k] = getInput<double>(mxGetField(prhs[3],k,"allN"));
   }

   arr(double) u = getArrayInput<double>(prhs[4]);
   arr(double) v = getArrayInput<double>(prhs[5]);

   plhs[0] = mxCreateNumericMatrix(xdim,ydim, mxDOUBLE_CLASS, mxREAL);
   plhs[1] = mxCreateNumericMatrix(xdim,ydim, mxDOUBLE_CLASS, mxREAL);
   arr(double) utotal = getArrayInput<double>(plhs[0]);
   arr(double) vtotal = getArrayInput<double>(plhs[1]);

   for (int i=0; i<xdim*ydim; i++)
   {
      utotal[i] = 10e10;
      vtotal[i] = 10e10;
   }

   // calculate the gradient information
   for (int k=0; k<K; k++)
   {
      arr(double) indices = indices_all[k];
      int N = indices_N[k];
      double uk = u[k];
      double vk = v[k];
      for (int i=0; i<N; i++)
      {
         int index = indices[i];
         utotal[index] = uk;
         vtotal[index] = vk;
      }
   }

   deallocate_memory(indices_all);
   deallocate_memory(indices_N);
}
