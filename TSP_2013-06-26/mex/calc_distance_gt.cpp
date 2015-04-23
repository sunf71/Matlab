// =============================================================================
// == calc_distance_gt.cp
// == --------------------------------------------------------------------------
// == A MEX interface to calculate the Boundary Recall Distance.
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

   checkInput(prhs, 0, MATRIX ,  "logical"); //gt_boundary
   checkInput(prhs, 1, MATRIX ,  "logical"); //sp_boundary
   checkInput(prhs, 2, VECTOR ,  "double"); //distances
   checkInput(prhs, 3, VECTOR ,  "int32"); //indicesx
   checkInput(prhs, 4, VECTOR ,  "int32"); //indicesy

   arr(bool) gt_boundary = getArrayInput<bool>(prhs[0]);
   int X = mxGetM(prhs[0]);
   int Y = mxGetN(prhs[0]);
   int N = X*Y;
   arr(bool) sp_boundary = getArrayInput<bool>(prhs[1]);
   arr(double) distances = getArrayInput<double>(prhs[2]);
   arr(int) indicesx = getArrayInput<int>(prhs[3]);
   arr(int) indicesy = getArrayInput<int>(prhs[4]);
   int D = mxGetNumberOfElements(prhs[2]);

   double total = 0;
   for (int i=0; i<N; i++)
   {
      if (gt_boundary[i])
      {
         int x = i%X;
         int y = i/X;
         for (int d=0; d<D; d++)
         {
            int nx = x + indicesx[d];
            int ny = y + indicesy[d];
            int ni = nx+ny*X;
            if (nx>=0 && nx<X && ny>=0 && ny<Y && sp_boundary[ni])
            {
               total += distances[d];
               break;
            }
         }
      }
   }

   plhs[0] = mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
   arr(double) temp = getArrayInput<double>(plhs[0]);
   temp[0] = total;
}
