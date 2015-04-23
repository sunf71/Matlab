// =============================================================================
// == SP_prop_init.cpp
// == --------------------------------------------------------------------------
// == A MEX interface to propagate labels according to some flow.
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
#include "linkedList.cpp"
#include "utils.h"

#define NUMARGS 7
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
   checkInput(prhs, 1, MATRIX ,  "int32" ); //label
   checkInput(prhs, 2, VECTOR ,  "double"); //meanx
   checkInput(prhs, 3, VECTOR ,  "double"); //meany
   checkInput(prhs, 4, VECTOR ,  "double"); //vx
   checkInput(prhs, 5, VECTOR ,  "double"); //vy
   checkInput(prhs, 6, MATRIX , "logical"); //mask

   int K = getInput<double>(prhs[0]);
   arr(int) label = getArrayInput<int>(prhs[1]);
   int X = mxGetM(prhs[1]);
   int Y = mxGetN(prhs[1]);
   int N = X*Y;
   double N2 = N*N;
   arr(double) meanx = getArrayInput<double>(prhs[2]);
   arr(double) meany = getArrayInput<double>(prhs[3]);
   arr(double) vx = getArrayInput<double>(prhs[4]);
   arr(double) vy = getArrayInput<double>(prhs[5]);
   arr(bool) mask = getArrayInput<bool>(prhs[6]);

   plhs[0] = mxCreateNumericMatrix(X,Y,mxINT32_CLASS,mxREAL);
   arr(int) newlabel = getArrayInput<int>(plhs[0]);

   memset(newlabel, -1, sizeof(int)*N);
   for (int i=0; i<N; i++)
      if (mask[i])
         newlabel[i] = K;

   arr( linkedList<int> ) possible_labels = allocate_memory< linkedList<int> >(N);

   for (int x=0; x<X; x++) for (int y=0; y<Y; y++)
   {
      int index = x + y*X;
      if (label[index]>=0)
      {
         int k = label[index];
         int newx = x + vx[k] + 0.5;
         int newy = y + vy[k] + 0.5;

         if (newx>=0 && newx<X && newy>=0 && newy<Y)
         {
            int newindex = newx + newy*X;
            possible_labels[newindex].addNodeEnd(k);
         }
      }
   }

   bool newRegion = false;
   for (int x=0; x<X; x++) for (int y=0; y<Y; y++)
   {
      int index = x + y*X;
      int count = possible_labels[index].getLength();
      if (count==1)
         newlabel[index] = possible_labels[index].popFirst();
      else if (count>1)
      {
         int closest_k = -1;
         double closest_distance = N2;
         while (!possible_labels[index].isempty())
         {
            int k = possible_labels[index].popFirst();
            double distance = pow(meanx[k] + vx[k] - x,2) + pow(meany[k] + vy[k] - y,2);
            if (distance<closest_distance || closest_k<0)
            {
               closest_distance = distance;
               closest_k = k;
            }
         }
         newlabel[index] = closest_k;
      }
      else if (mask[index])
         newRegion = true;
   }
   if (newRegion)
      K++;

   // labels propagated, now enforce connectivity
   arr(bool) done = allocate_memory<bool>(N, false);
   arr( linkedList< linkedList<int>* >) indices = allocate_memory< linkedList< linkedList<int>* > >(K);
   linkedList<int> explore;
   for (int x=0; x<X; x++) for (int y=0; y<Y; y++)
   {
      int index = x + y*X;
      if (!done[index] && newlabel[index]>=0)
      {
         done[index] = true;
         int curLabel = newlabel[index];
         explore.addNodeEnd(index);
         linkedList<int>* thisList = new linkedList<int>();
         indices[curLabel].addNodeEnd(thisList);

         while (!explore.isempty())
         {
            index = explore.popFirst();
            thisList->addNodeEnd(index);
            int xi = index%X;
            int yi = index/X;
            if (xi>0 && !done[index-1] && newlabel[index-1]==curLabel)
            {
               explore.addNodeEnd(index-1);
               done[index-1] = true;
            }
            if (yi>0 && !done[index-X] && newlabel[index-X]==curLabel)
            {
               explore.addNodeEnd(index-X);
               done[index-X] = true;
            }
            if (xi<X-1 && !done[index+1] && newlabel[index+1]==curLabel)
            {
               explore.addNodeEnd(index+1);
               done[index+1] = true;
            }
            if (yi<Y-1 && !done[index+X] && newlabel[index+X]==curLabel)
            {
               explore.addNodeEnd(index+X);
               done[index+X] = true;
            }
         }
      }
      else if (!done[index] && newlabel[index]<0)
         done[index] = true;
   }

   int curK = K;
   for (int k=0; k<curK; k++)
   {
      if (indices[k].getLength()>1)
      {
         // find maximum length one first
         linkedListNode< linkedList<int>* >* theListNode = indices[k].getFirst();
         int maxLength = 0;
         int maxLengthi = -1;
         int i = 0;
         while (theListNode!=NULL)
         {
            int theLength = theListNode->getData()->getLength();
            if (theLength > maxLength)
            {
               maxLength = theLength;
               maxLengthi = i;
            }
            theListNode = theListNode->getNext();
            i++;
         }

         // mark all the other ones as not being done
         i = 0;
         theListNode = indices[k].getFirst();
         while (theListNode!=NULL)
         {
            if (i!=maxLengthi)
            {
               linkedList<int>* theList = theListNode->getData();
               linkedListNode<int>* theNode = theList->getFirst();

               if (theList->getLength()<20)
               {
                  while (theNode!=NULL)
                  {
                     int index = theNode->getData();
                     done[index] = false;
                     theNode = theNode->getNext();
                  }
               }
               else
               {
                  while (theNode!=NULL)
                  {
                     int index = theNode->getData();
                     done[index] = true;
                     theNode = theNode->getNext();
                     newlabel[index] = K;
                  }
                  K++;
               }
            }
            theListNode = theListNode->getNext();
            i++;
         }
      }
      if (indices[k].getLength()>0)
      {
         linkedListNode< linkedList<int>* >* theListNode = indices[k].getFirst();
         while (theListNode!=NULL)
         {
            delete theListNode->getData();
            theListNode = theListNode->getNext();
         }
      }
   }


   bool any_notdone = true;
   int count = 0;
   while (any_notdone)
   {
      any_notdone = false;
      // we can either just set to a neighboring K, or create new ones
      // this sets to a neighboring k
      for (int x=0; x<X; x++) for (int y=0; y<Y; y++)
      {
         int index = x + y*X;
         if (!done[index])
         {
            done[index] = true;
            if (x<X-1 && done[index+1] && newlabel[index+1]>=0)
               newlabel[index] = newlabel[index+1];
            else if (y<Y-1 && done[index+X] && newlabel[index+X]>=0)
               newlabel[index] = newlabel[index+X];
            else if (x>0 && done[index-1] && newlabel[index-1]>=0)
               newlabel[index] = newlabel[index-1];
            else if (y>0 && done[index-X] && newlabel[index-X]>=0)
               newlabel[index] = newlabel[index-X];
            else
               done[index] = false;
         }
      }
      for (int x=X-1; x>=0; x--) for (int y=Y-1; y>=0; y--)
      {
         int index = x + y*X;
         if (!done[index])
         {
            done[index] = true;
            if (x<X-1 && done[index+1] && newlabel[index+1]>=0)
               newlabel[index] = newlabel[index+1];
            else if (y<Y-1 && done[index+X] && newlabel[index+X]>=0)
               newlabel[index] = newlabel[index+X];
            else if (x>0 && done[index-1] && newlabel[index-1]>=0)
               newlabel[index] = newlabel[index-1];
            else if (y>0 && done[index-X] && newlabel[index-X]>=0)
               newlabel[index] = newlabel[index-X];
            else
            {
               done[index] = false;
               any_notdone = true;
            }
         }
      }
      count++;
      if (count>10)
      {
         break;
         //mexErrMsgTxt("SP_prop_init exceeded tries\n");
      }
   }

   for (int i=0; i<N; i++) if (!done[i])
   {
      if (mask[i])
         mexErrMsgTxt("SP_prop_init exceeded tries\n");
      else
         newlabel[i] = -1;
   }

   deallocate_memory(done);
   deallocate_memory(indices);

   deallocate_memory(possible_labels);
}
