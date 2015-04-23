#include "mex.h"
#include "helperMEX.h"
#include "is_border.cpp"

#define NUMARGS 1
#define NUMOUT 1
//void extend_image(double xdimDouble, double ydimDouble, double* imageIn,
//   double* maskIn, bool linearExtension, double* distance,
//   double* closest_curve_x, double* closest_curve_y, double* farthest_size,
//   double* farthest_x, double* farthest_y, double* imageOut)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *inptr1;
  bool *outptr1;
  bool inparam5;
  int mrows, ncols;

  // Check for proper number of arguments
  if (nrhs != NUMARGS) {
      mexErrMsgTxt("Incorrect number of input arguments required.");
  } else if (nlhs > NUMOUT) {
      mexErrMsgTxt("Too many output arguments expected.");
  }

  checkInput(prhs, 0, 'm');

  //Create matrix for the return argument
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  plhs[0] = mxCreateLogicalMatrix(mrows,ncols);
  //Assign pointers to input and output
  inptr1 = (double *) mxGetData(prhs[0]);

  outptr1 = (bool *) mxGetData(plhs[0]);
  /* Retrieve the actual values of the
     input arguments. */

  /* Call the subroutine. */
  // SUBROUTINE() is defined in the
  //FunctionNameIMPORT.cpp
  // or FunctionNameEXPORT.cpp file
  is_border_vals(inptr1, mrows, ncols, outptr1);
}
