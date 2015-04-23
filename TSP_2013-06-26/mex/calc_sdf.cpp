// =============================================================================
// == calc_sdf.cpp
// == --------------------------------------------------------------------------
// == The MEX interface file to calculate a signed distance function
// == --------------------------------------------------------------------------
// == Copyright 2011. MIT. All Rights Reserved.
// == Written by Jason Chang 06-13-2011
// =============================================================================

#include "binTree.cpp"
#include "distancePixel.cpp"
#include "math.h"
#include "linkedList.cpp"


// --------------------------------------------------------------------------
// -- calc_dist
// --   calculates distance to zero level set at far point (x,y)
// --
// --   parameters:
// --     - (x,y) : the coordinate where the distance is needed
// --     - (xNew, yNew) : the newly accepted distance that touches (x,y)
// --     - unmarked : the tree of points that have a finite distance but
// --         have not been accepted
// --     - accepted : the 2D boolean matrix saying if the pixel is accepted
// --     - tree_ptr : a 2D matrix of pointers pointing to the node in the
// --         unmarked tree (if exists) or NULL (if doesn't exist)
// --     - band_size : the narrow band size
// --
// --   return parameters:
// --     - distance : the SDF, updated (x,y)
// --------------------------------------------------------------------------
inline void calc_dist(const int x, const int y, const int xNew, const int yNew,
   int xdim, int ydim, binTree<distancePixel> &unmarked, bool* accepted,
   binNode<distancePixel>** tree_ptr, double band_size, double* distance)
{
   // do x direction
   double a = 0;
   double b = 0;
   double c = -1;
   
   bool checkbck = (x-1>=0) && accepted[(x-1)*ydim+y];
   bool checkfwd = (x+1<xdim) && accepted[(x+1)*ydim+y];
   double dist = 10*band_size;
   if (checkbck && !checkfwd)
      dist = distance[(x-1)*ydim+y];
   else if (!checkbck && checkfwd)
      dist = distance[(x+1)*ydim+y];
   else if (checkbck && checkfwd)
      dist = (distance[(x+1)*ydim+y] > distance[(x-1)*ydim+y]) ? distance[(x-1)*ydim+y] : distance[(x+1)*ydim+y];
   if (dist != 10*band_size)
   {
      a += 1;
      b += -2*dist;
      c += dist*dist;
   }

   // do y direction
   checkbck = (y-1>=0) && accepted[x*ydim+y-1];
   checkfwd = (y+1<ydim) && accepted[x*ydim+y+1];
   dist = 10*band_size;
   if (checkbck && !checkfwd)
      dist = distance[x*ydim+y-1];
   else if (!checkbck && checkfwd)
      dist = distance[x*ydim+y+1];
   else if (checkbck && checkfwd)
      dist = (distance[x*ydim+y+1] > distance[x*ydim+y-1]) ? distance[x*ydim+y-1] : distance[x*ydim+y+1];
   if (dist != 10*band_size)
   {
      a += 1;
      b += -2*dist;
      c += dist*dist;
   }

   double sign = (distance[x*ydim+y] > 0) ? 1 : -1;
   dist = (-b + sign * sqrt(b*b - 4*a*c)) / (2*a);
   //if (dist < 0) dist = -dist;

   if (fabs(dist) < fabs(distance[x*ydim+y]) && fabs(dist)<=band_size)
   {
      // add it to the tree because it had never been added before
      if (tree_ptr[x*ydim+y] == NULL)
         tree_ptr[x*ydim+y] = unmarked.addLeaf(distancePixel(dist, x, y));
      // modify it from the tree
      else
         unmarked.modifyLeaf(tree_ptr[x*ydim+y], distancePixel(dist, x, y));
      distance[x*ydim+y] = dist;
   }
}

// --------------------------------------------------------------------------
// -- calc_sdf
// --   calculates a signed distance function for distanceIn (which only must
// -- contain the correct sign). Pixel i is only updated if maskIn[i]=true,
// -- and in a narrow band region if size band_size.
// --
// --   parameters:
// --     - (xdimDouble, ydimDouble) : the dimensions of the image
// --     - distanceIn : the new level set function to calculate the SDF of.
// --         Only the signs of this matter.
// --     - maskIn : the boolean array indicating which pixels to update
// --     - oldDistance : the old values of the SDF
// --     - band_size : the size of the narrow band
// --
// --   return parameters:
// --     - distance : the calculated SDF
// --------------------------------------------------------------------------
inline void calc_sdf(double xdimDouble, double ydimDouble, double* distanceIn,
   bool* maskIn, double* oldDistance, double band_size,
   double* &distance)
{
   int xdim = (int)xdimDouble;
   int ydim = (int)ydimDouble;

   binTree<distancePixel> unmarked;

   bool* accepted = new bool[xdim*ydim];
   binNode<distancePixel>** tree_ptr = new binNode<distancePixel>*[xdim*ydim];
   int index;

   for (int x=0; x<xdim; x++) for (int y=0; y<ydim; y++)
   {
      index = x*ydim+y;
      tree_ptr[index] = NULL;
      if (!maskIn[index])
      {
         accepted[index] = true;
         distance[index] = oldDistance[index];
      }
      else if (((distanceIn[index]>=0) && ((x>0 && distanceIn[index-ydim]<0) ||
                                           (x<xdim-1 && distanceIn[index+ydim]<0) ||
                                           (y>0 && distanceIn[index-1]<0) ||
                                           (y<ydim-1 && distanceIn[index+1]<0))) ||
               ((distanceIn[index]<0 ) && ((x>0 && distanceIn[index-ydim]>=0) ||
                                           (x<xdim-1 && distanceIn[index+ydim]>=0) ||
                                           (y>0 && distanceIn[index-1]>=0) ||
                                           (y<ydim-1 && distanceIn[index+1]>=0))))
      {
         accepted[index] = true;
         distance[index] = 0;

         tree_ptr[index] = unmarked.addLeaf( distancePixel(0.5 * (distanceIn[index]>=0 ? 1 : -1), x, y) );
      }
      else
      {
         if (distanceIn[index]<0)
            distance[index] = -band_size;
         else
            distance[index] = band_size;
         accepted[index] = false;
      }
   }

   // recursively pick off the smallest distance until all pixels are done
   while (unmarked.treeNotEmpty())
   {
      distancePixel new_point = unmarked.pickOffRoot();
      int x = new_point.coords.xpix;
      int y = new_point.coords.ypix;
      accepted[x*ydim+y] = true;
      tree_ptr[x*ydim+y] = NULL;
      distance[x*ydim+y] = new_point.distance;

      if (x+1<xdim && !accepted[(x+1)*ydim+y] ) calc_dist(x+1,y,x,y, xdim, ydim, unmarked, accepted, tree_ptr, band_size, distance);
      if (x-1>=0   && !accepted[(x-1)*ydim+y] ) calc_dist(x-1,y,x,y, xdim, ydim, unmarked, accepted, tree_ptr, band_size, distance);
      if (y+1<ydim && !accepted[x*ydim+(y+1)] ) calc_dist(x,y+1,x,y, xdim, ydim, unmarked, accepted, tree_ptr, band_size, distance);
      if (y-1>=0   && !accepted[x*ydim+(y-1)] ) calc_dist(x,y-1,x,y, xdim, ydim, unmarked, accepted, tree_ptr, band_size, distance);
   }

   delete[] tree_ptr;
   delete[] accepted;
}
