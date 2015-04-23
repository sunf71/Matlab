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
#include "utils.h"

#define NUMARGS 2
#define NUMOUT 2

inline bool is_border(int x, int y, int X, int Y, arr(int) label)
{
   int index = x + y*X;
   bool border = false;
   int cur_label = label[index];
   if (x>0) border = border || (cur_label!=label[index-1]);
   if (y>0) border = border || (cur_label!=label[index-X]);
   if (x<X-1) border = border || (cur_label!=label[index+1]);
   if (y<Y-1) border = border || (cur_label!=label[index+X]);
   return border;
}

inline bool is_border_f(int x, int y, int X, int Y, arr(int) label, int& border_label)
{
   int index = x + y*X;
   bool border = false;
   int cur_label = label[index];
   if (x>=X-1)
      border_label = -1;
   else
      border_label = label[index+1];
   return (x>=X-1) || (cur_label!=label[index+1]);
}
inline bool is_border_b(int x, int y, int X, int Y, arr(int) label, int& border_label)
{
   int index = x + y*X;
   bool border = false;
   int cur_label = label[index];
   if (x<=0)
      border_label = -1;
   else
      border_label = label[index-1];
   return (x<=0) || (cur_label!=label[index-1]);
}
inline bool is_border_d(int x, int y, int X, int Y, arr(int) label, int& border_label)
{
   int index = x + y*X;
   bool border = false;
   int cur_label = label[index];
   if (y>=Y-1)
      border_label = -1;
   else
      border_label = label[index+X];
   return (y>=Y-1) || (cur_label!=label[index+X]);
}
inline bool is_border_u(int x, int y, int X, int Y, arr(int) label, int& border_label)
{
   int index = x + y*X;
   bool border = false;
   int cur_label = label[index];
   if (y<=0)
      border_label = -1;
   else
      border_label = label[index-X];
   return (y<=0) || (cur_label!=label[index-X]);
}

inline void add_neighbor(int k, int border_label, int K, arr(bool) neighbors, arr(int) neighbors_N)
{
   if (border_label>=0)
   {
      if (!neighbors[k+border_label*K])
      {
         neighbors[k+border_label*K] = true;
         neighbors_N[k]++;
      }
   }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   // Check for proper number of arguments
   if (nrhs != NUMARGS) {
         mexErrMsgTxt("Incorrect number of input arguments required.");
   } else if (nlhs > NUMOUT) {
         mexErrMsgTxt("Too many output arguments expected.");
   }

   checkInput(prhs, 0, SCALAR ,  "double"); //K
   checkInput(prhs, 1, MATRIX ,  "int32"); //label
mexPrintf("a\n");
   int K = getInput<double>(prhs[0]);
   arr(int) label = getArrayInput<int>(prhs[1]);
   int X = mxGetM(prhs[1]);
   int Y = mxGetN(prhs[1]);
mexPrintf("b\n");
   // outputs
   arr( bool ) neighbors = allocate_memory<bool>(K*K,0);
   arr( int ) neighbors_N = allocate_memory<int>(K,0);
   arr( int ) indices_N = allocate_memory<int>(K);
   arr( int ) indices_allN = allocate_memory<int>(K);
   arr( linkedList<int> ) indices_all = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_a = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_c = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_f = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_b = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_d = allocate_memory< linkedList<int> >(K);
   arr( linkedList<int> ) indices_u = allocate_memory< linkedList<int> >(K);
mexPrintf("c\n");
   for (int x=0; x<X; x++)
   {
      for (int y=0; y<Y; y++)
      {
         int index = x + y*X;
         int k = label[index];
         if (k>=0)
         {
            int border_label = -1;

            indices_all[k].addNode(index);
            if (!is_border(x,y,X,Y,label))
               indices_a[k].addNode(index);
            else
               indices_c[k].addNode(index);

            if (is_border_f(x,y,X,Y,label, border_label))
            {
               indices_f[k].addNode(index);
               add_neighbor(k, border_label, K, neighbors, neighbors_N);
            }
            if (is_border_b(x,y,X,Y,label, border_label))
            {
               indices_b[k].addNode(index);
               add_neighbor(k, border_label, K, neighbors, neighbors_N);
            }
            if (is_border_d(x,y,X,Y,label, border_label))
            {
               indices_d[k].addNode(index);
               add_neighbor(k, border_label, K, neighbors, neighbors_N);
            }
            if (is_border_u(x,y,X,Y,label, border_label))
            {
               indices_u[k].addNode(index);
               add_neighbor(k, border_label, K, neighbors, neighbors_N);
            }

            indices_N[k]++;
            indices_allN[k]++;
         }
      }
   }
mexPrintf("d\n");
   // indices
   const char* indices_names[9] = {"N", "allN", "all", "a", "c", "f", "b", "d", "u"};
   plhs[0] = mxCreateStructMatrix(K,1,9,indices_names);
   for (int k=0; k<K; k++)
   {
      mxWriteField(plhs[0],k,"N",   mxWriteScalar((double)(indices_N[k])));
      mxWriteField(plhs[0],k,"allN",mxWriteScalar((double)(indices_allN[k])));
      mxWriteField(plhs[0],k,"all", mxWriteVector(1,indices_all[k].getLength(),indices_all[k].getArray()));
      mxWriteField(plhs[0],k,"a",   mxWriteVector(1,indices_a[k].getLength(),indices_a[k].getArray()));
      mxWriteField(plhs[0],k,"c",   mxWriteVector(1,indices_c[k].getLength(),indices_c[k].getArray()));
      mxWriteField(plhs[0],k,"f",   mxWriteVector(1,indices_f[k].getLength(),indices_f[k].getArray()));
      mxWriteField(plhs[0],k,"b",   mxWriteVector(1,indices_b[k].getLength(),indices_b[k].getArray()));
      mxWriteField(plhs[0],k,"d",   mxWriteVector(1,indices_d[k].getLength(),indices_d[k].getArray()));
      mxWriteField(plhs[0],k,"u",   mxWriteVector(1,indices_u[k].getLength(),indices_u[k].getArray()));
   }
mexPrintf("e\n");
   // neighbors
   plhs[1] = mxCreateCellMatrix(K,1);
   for (int k=0; k<K; k++)
   {
      mxSetCell(plhs[1],k, mxCreateNumericMatrix(1,neighbors_N[k], mxDOUBLE_CLASS, mxREAL));
      arr(double) neighbors_k = getArrayInput<double>(mxGetCell(plhs[1],k));
      int count = 0;
      for (int k2=0; k2<K; k2++)
         if (neighbors[k+k2*K])
            neighbors_k[count++] = k2;
   }
mexPrintf("f\n");
   deallocate_memory(neighbors);
   deallocate_memory(neighbors_N);
   deallocate_memory(indices_N);
   deallocate_memory(indices_allN);
   deallocate_memory(indices_all);
   deallocate_memory(indices_a);
   deallocate_memory(indices_c);
   deallocate_memory(indices_f);
   deallocate_memory(indices_b);
   deallocate_memory(indices_d);
   deallocate_memory(indices_u);
}
