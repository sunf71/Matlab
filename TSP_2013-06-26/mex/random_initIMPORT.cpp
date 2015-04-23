// =============================================================================
// == calc_sdfIMPORT.cpp
// == --------------------------------------------------------------------------
// == The MEX interface file to calculate a signed distance function
// == --------------------------------------------------------------------------
// == Copyright 2011. MIT. All Rights Reserved.
// == Written by Jason Chang 06-13-2011
// =============================================================================

#include "mex.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <time.h>
#include "debugMEX.h"
#include "helperMEX.h"
#include "matrix.h"
#include "linkedList.cpp"

#define NUMARGS 3
#define NUMOUT 1


// find_valid_noiseIMPORT(phi, mask, cc, logpim_p, logpim_m, Eright_same, Eright_diff, Edown_same, Edown_diff, FG_topologies, BG_topologies, allowInsert, allowRemove)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, SCALAR,  "double"); //xdim
   checkInput(prhs, 1, SCALAR,  "double"); //ydim
   checkInput(prhs, 2, VECTOR,  "double"); //centers

   int xdim = getInput<double>(prhs[0]);
   int ydim = getInput<double>(prhs[1]);
   int N = xdim*ydim;
   arr(double) centers = getArrayInput<double>(prhs[2]);
   int K = mxGetNumberOfElements(prhs[2]);

   plhs[0] = mxCreateNumericMatrix(xdim, ydim, mxDOUBLE_CLASS, mxREAL);
   arr(double) z = getArrayInput<double>(plhs[0]);
   for (int i=0; i<N; i++)
      z[i] = -1;

   linkedList<unsigned int> indices;
   linkedList<unsigned int> indices_z;
   for (int k=0; k<K; k++)
   {
      indices.addNodeEnd(centers[k]-1);
      indices_z.addNodeEnd(k+1);
      z[(int)(centers[k])-1] = k+1;
   }

   while (!(indices.isempty()))
   {
      unsigned int index = indices.popFirst();
      unsigned int index_z = indices_z.popFirst();

      // add neighbors
      int x = index%xdim;
      int y = index/xdim;
      if (x>0 && z[index-1]<0)
      {
         z[index-1] = index_z;
         indices.addNodeEnd(index-1);
         indices_z.addNodeEnd(index_z);
      }
      if (y>0 && z[index-xdim]<0)
      {
         z[index-xdim] = index_z;
         indices.addNodeEnd(index-xdim);
         indices_z.addNodeEnd(index_z);
      }
      if (x<xdim-1 && z[index+1]<0)
      {
         z[index+1] = index_z;
         indices.addNodeEnd(index+1);
         indices_z.addNodeEnd(index_z);
      }
      if (y<ydim-1 && z[index+xdim]<0)
      {
         z[index+xdim] = index_z;
         indices.addNodeEnd(index+xdim);
         indices_z.addNodeEnd(index_z);
      }
   }
}
