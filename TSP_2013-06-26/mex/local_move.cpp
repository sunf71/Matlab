// =============================================================================
// == local_move.cpp
// == --------------------------------------------------------------------------
// == A MEX interface to perform local TSP moves with flow estimation.
// == See m files for calling convention.
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================


#include "mex.h"
#include "IMG.h"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )
{
   if (nrhs!=2) mexErrMsgTxt("Only the IMG structure and # of its expected");
   if (mxGetNumberOfElements(prhs[0])!=1) mexErrMsgTxt("One structure expected.");

   IMG sp_img;
   sp_img.mxReadIMG(prhs[0]);

   int its = getInput<double>(prhs[1]);
   bool converged = false;
   for (int i=0; i<its; i++)
   {
      if (i%20==0 || i==its-1)
      {
         sp_img.Flow_QP2();
         converged = true;
         if (!sp_img.move_local_IMG())
            break;
      }
      else
      {
         if (!converged)
            converged = sp_img.move_local_IMG();
      }
   }
   sp_img.mxWriteIMG(plhs,prhs[0]);
}
