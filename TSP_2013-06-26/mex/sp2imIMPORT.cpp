// =============================================================================
// == sp2imIMPORT.cpp
// == --------------------------------------------------------------------------
// == A MEX interface to convert SP data to an image
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

#define NUMARGS 5
#define NUMOUT 1


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
   checkInput(prhs, 3, CELL              ); //indices
   checkInput(prhs, 4, VECTOR ,  "double"); //data


   int K = getInput<double>(prhs[0]);
   int xdim = getInput<double>(prhs[1]);
   int ydim = getInput<double>(prhs[2]);

   arr( arr(int) ) indices_all = allocate_memory<arr(int)>(K);
   arr(int) indices_N = allocate_memory<int>(K);
   for (int k=0; k<K; k++)
   {
      if (mxGetCell(prhs[3],k))
      {
         if (mxGetNumberOfElements(mxGetCell(prhs[3],k))==1)
            checkInput(mxGetCell(prhs[3],k),SCALAR,"uint32");
         else
            checkInput(mxGetCell(prhs[3],k),VECTOR,"uint32");
         indices_all[k] = getArrayInput<int>(mxGetCell(prhs[3],k));
         indices_N[k] = mxGetNumberOfElements(mxGetCell(prhs[3],k));
      }
      else
      {
         indices_all[k] = NULL;
         indices_N[k] = 0;
      }
   }

   arr(double) spdata = getArrayInput<double>(prhs[4]);

   plhs[0] = mxCreateNumericMatrix(xdim,ydim, mxDOUBLE_CLASS, mxREAL);
   arr(double) imdata = getArrayInput<double>(plhs[0]);

   for (int i=0; i<xdim*ydim; i++)
   {
      imdata[i] = 10e10;
   }

   // calculate the gradient information
   for (int k=0; k<K; k++)
   {
      arr(int) indices = indices_all[k];
      int N = indices_N[k];
      double spdatai = spdata[k];
      for (int i=0; i<N; i++)
      {
         int index = indices[i];
         imdata[index] = spdatai;
      }
   }

   deallocate_memory(indices_all);
   deallocate_memory(indices_N);
}
