// =============================================================================
// == calc_sdfIMPORT.cpp
// == --------------------------------------------------------------------------
// == The MEX interface file to calculate a signed distance function
// == --------------------------------------------------------------------------
// == Copyright 2011. MIT. All Rights Reserved.
// == Written by Jason Chang 06-13-2011
// =============================================================================

#include "mex.h"
#include "calc_sdf.cpp"
#include <iostream>
#include "helperMEX.h"

#define NUMARGS 4
#define NUMOUT 1
// sdf = calc_sdfIMPORT(regions, mask, oldsdf, band_size);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   int mrows, ncols;

   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, MATRIX ,  "int32"); //label

   //Create matrix for the return argument
   mrows = mxGetM(prhs[2]);
   ncols = mxGetN(prhs[2]);
   plhs[0] = mxCreateNumericMatrix(mrows,ncols, mxDOUBLE_CLASS, mxREAL);

   //Assign pointers to input and output
   arr(int) regions = getArrayInput<int>(mxGetData(prhs[0]));
   arr(double) sdf = getArrayInput<double>(mxGetData(plhs[0]));

   calc_sdf(ncols, mrows, regions, sdf);
}
