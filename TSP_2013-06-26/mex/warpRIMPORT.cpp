// =============================================================================
// == warpRIMPORT.cpp
// == --------------------------------------------------------------------------
// == A MEX interface to reverse warp an image according to flow
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
#include "debugMEX.h"
#include "helperMEX.h"


#define NUMARGS 3
#define NUMOUT 1

void calc_weights(double x, double y,
                  double &UL, double &UR, double &DL, double &DR)
{
   x = x - (int)(x);
   y = y - (int)(y);
   if (x<0) x+=1;
   if (y<0) y+=1;

   UL = 1-x;
   DL = 1-x;
   UR = x;
   DR = x;
   UL *= (1-y);
   UR *= (1-y);
   DL *= y;
   DR *= y;
}



// find_valid_noiseIMPORT(phi, mask, cc, logpim_p, logpim_m, Eright_same, Eright_diff, Edown_same, Edown_diff, FG_topologies, BG_topologies, allowInsert, allowRemove)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs < NUMARGS ) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   }/* else if (nlhs) {
         mexErrMsgTxt("Number of output warpings must be same as input.");
   }*/

   int Nimages = nlhs;

   checkInput(prhs[0], MATRIX, "double"); //fx
   checkInput(prhs[1], MATRIX, "double"); //fy

   int xdim = mxGetM(prhs[0]);
   int ydim = mxGetN(prhs[0]);


   checkInput(prhs[2], MATRIX, "double");

   int N = xdim*ydim;
   arr(double) fy = getArrayInput<double>(prhs[0]);
   arr(double) fx = getArrayInput<double>(prhs[1]);

   arr(double) input = getArrayInput<double>(prhs[2]);
   const int *dims = mxGetDimensions(prhs[2]);
   int ndim;
   if (mxGetNumberOfDimensions(prhs[2])==2)
      ndim = 1;
   else
      ndim = dims[2];
   plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[2]), dims, mxDOUBLE_CLASS, mxREAL);
   plhs[1] = mxCreateNumericMatrix(xdim,ydim, mxDOUBLE_CLASS, mxREAL);
   arr(double) output = getArrayInput<double>(plhs[0]);
   arr(double) total_weight = getArrayInput<double>(plhs[1]);

   memset(total_weight, 0, sizeof(double)*N);
   //arr(double) total_weight = allocate_memory<double>(N, 0);

   double UL, UR, DL, DR, UL2, UR2, DL2, DR2;
   for (int i=0; i<N; i++)
   {
      int x = i%xdim;
      int y = i/xdim;

      double nx = fx[i] + x;
      double ny = fy[i] + y;

      int ni = (int)(nx)+(int)(ny)*xdim;

      calc_weights(nx, ny, UL, UR, DL, DR);

      if (UL<0 || UR<0 || DL<0 || DR<0)
         mexErrMsgTxt("WTF!\n");

      // UL
      bool UL_good = (nx>=0 && nx<xdim && ny>=0 && ny<ydim);
      bool UR_good = (nx>=-1 && nx<xdim-1 && ny>=0 && ny<ydim);
      bool DL_good = (nx>=0 && nx<xdim && ny>=-1 && ny<ydim-1);
      bool DR_good = (nx>=-1 && nx<xdim-1 && ny>=-1 && ny<ydim-1);

      UL2 = UL*UL;
      UR2 = UR*UR;
      DL2 = DL*DL;
      DR2 = DR*DR;

      if (UL_good)
      {
         total_weight[ni] += UL;
//         total_weight2[ni] += UL2;
      }
      if (UR_good)
      {
         total_weight[ni+1] += UR;
//         total_weight2[ni+1] += UR2;
      }
      if (DL_good)
      {
         total_weight[ni+xdim] += DL;
//         total_weight2[ni+xdim] += DL2;
      }
      if (DR_good)
      {
         total_weight[ni+xdim+1] += DR;
//         total_weight2[ni+xdim+1] += DR2;
      }

      for (int d=0, ind=i, nind=ni; d<ndim; d++, ind+=N, nind+=N)
      {
         double data = input[ind];
         if (UL_good) output[nind] += data * UL;
         if (UR_good) output[nind+1] += data * UR;
         if (DL_good) output[nind+xdim] += data * DL;
         if (DR_good) output[nind+xdim+1] += data * DR;
      }
   }

   for (int i=0; i<N; i++)
   {
      double this_total = total_weight[i];
      double this_totalsq = this_total*this_total;

      if (this_total>0)
      {
         for (int d=0, ind=i; d<ndim; d++, ind+=N)
            output[ind] /= this_total;
      }
   }
}
