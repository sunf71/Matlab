// =============================================================================
// == is_border.cpp
// == --------------------------------------------------------------------------
// == Some functions to help in declaring boundaries
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================

#include "math.h"

void is_border_vals(double* level, int xdim, int ydim, bool* border);
void is_border(double* level, int xdim, int ydim, bool* border);
void is_border_corner(double* level, int xdim, int ydim, bool* border);

void is_border_vals(double* level, int xdim, int ydim, bool* border)
{
	for (int x=0; x<xdim; x++) for (int y=0; y<ydim; y++)
	{
		int x1 = (x-1<0) ? 0 : x-1;
		int y1 = (y-1<0) ? 0 : y-1;
		int x2 = (x+1<xdim) ? x+1 : xdim-1;
		int y2 = (y+1<ydim) ? y+1 : ydim-1;

		double val = level[y*xdim+x];
		if (val != level[y*xdim+x1] || val != level[y1*xdim+x] ||
		    val != level[y*xdim+x2] || val != level[y2*xdim+x])
			border[y*xdim+x] = true;
		else
			border[y*xdim+x] = false;
	}
}

void is_border(double* level, int xdim, int ydim, bool* border)
{
   // handle the image corners first
   int i, offset;
   bool s, s1, s2, s3, s4;
   i = 0;
   s = level[i]>0;
   s1 = level[i+1]>0;
   s2 = level[i+xdim]>0;
   if (s!=s1 || s!=s2)
      border[0] = true;
   else
      border[0] = false;

   i = xdim-1;
   s = level[i]>0;
   s1 = level[i-1]>0;
   s2 = level[i+xdim]>0;
   if (s!=s1 || s!=s2)
      border[i] = true;
   else
      border[i] = false;

   i = xdim*(ydim-1);
   s = level[i]>0;
   s1 = level[i+1]>0;
   s2 = level[i-xdim]>0;
   if (s!=s1 || s!=s2)
      border[i] = true;
   else
      border[i] = false;

   i = xdim*ydim-1;
   s = level[i]>0;
   s1 = level[i-1]>0;
   s2 = level[i-xdim]>0;
   if (s!=s1 || s!=s2)
      border[i] = true;
   else
      border[i] = false;

   // handle the image borders first
   s = level[0]>0;
   s3 = level[1]>0;
   for (int x=1; x<xdim-1; x++)
   {
      s1 = s;
      s2 = level[x+xdim]>0;
      s = s3;
      s3 = level[x+1]>0;
      if (s!=s1 || s!=s2 || s!=s3)
         border[x] = true;
      else
         border[x] = false;
   }

   offset = xdim*(ydim-1);
   s = level[offset]>0;
   s3 = level[offset+1]>0;
   for (int x=1; x<xdim-1; x++)
   {
      s1 = s;
      s2 = level[offset+x-xdim]>0;
      s = s3;
      s3 = level[offset+x+1]>0;
      if (s!=s1 || s!=s2 || s!=s3)
         border[offset+x] = true;
      else
         border[offset+x] = false;
   }

   s = level[0]>0;
   s3 = level[xdim]>0;
   for (int y=1; y<ydim-1; y++)
   {
      s1 = s;
      s2 = level[y*xdim+1]>0;
      s = s3;
      s3 = level[y*xdim+xdim]>0;
      if (s!=s1 || s!=s2 || s!=s3)
         border[y*xdim] = true;
      else
         border[y*xdim] = false;
   }

   s = level[xdim-1]>0;
   s3 = level[2*xdim-1]>0;
   for (int y=1; y<ydim-1; y++)
   {
      s1 = s;
      s2 = level[y*xdim+xdim-2]>0;
      s = s3;
      s3 = level[y*xdim+xdim-1+xdim]>0;
      if (s!=s1 || s!=s2 || s!=s3)
         border[y*xdim+xdim-1] = true;
      else
         border[y*xdim+xdim-1] = false;
   }

   int newoffset;
   for (int y=1; y<ydim-1; y++)
   {
      offset = y*xdim;

      s1 = level[offset-1]>0;
      s = level[offset]>0;
      s3 = level[offset+1]>0;
      for (int x=1; x<xdim-1; x++)
      {
         s1 = s;
         s = s3;

         newoffset = offset + x;
         s2 = level[newoffset+xdim]>0;
         s3 = level[newoffset+1]>0;
         s4 = level[newoffset-xdim]>0;

         if (s!=s1 || s!=s2 || s!=s3 || s!=s4)
            border[newoffset] = true;
         else
            border[newoffset] = false;
      }
   }
}





void is_border_corner(double* level, int xdim, int ydim, bool* border)
{
   // handle the image corners first
   int i, offset;
   bool s, s1, s2, s3, s4, s5, s6, s7, s8;
   i = 0;
   s = level[i]>0;
   s1 = level[i+1]>0;
   s2 = level[i+xdim]>0;
   s3 = level[i+xdim+1]>0;
   if (s!=s1 || s!=s2 || s!=s3)
      border[0] = true;
   else
      border[0] = false;

   i = xdim-1;
   s = level[i]>0;
   s1 = level[i-1]>0;
   s2 = level[i+xdim]>0;
   s3 = level[i+xdim-1]>0;
   if (s!=s1 || s!=s2 || s!=s3)
      border[i] = true;
   else
      border[i] = false;

   i = xdim*(ydim-1);
   s = level[i]>0;
   s1 = level[i+1]>0;
   s2 = level[i-xdim]>0;
   s3 = level[i-xdim+1]>0;
   if (s!=s1 || s!=s2 || s!=s3)
      border[i] = true;
   else
      border[i] = false;

   i = xdim*ydim-1;
   s = level[i]>0;
   s1 = level[i-1]>0;
   s2 = level[i-xdim]>0;
   s3 = level[i-xdim-1]>0;
   if (s!=s1 || s!=s2 || s!=s3)
      border[i] = true;
   else
      border[i] = false;

   // handle the image borders first
   s = level[0]>0;
   s3 = level[xdim]>0;
   s4 = level[xdim+1]>0;
   s5 = level[1]>0;
   for (int x=1; x<xdim-1; x++)
   {
      s1 = s;
      s2 = s3;
      s = s5;
      s3 = s4;
      s4 = level[x+xdim+1]>0;
      s5 = level[x+1]>0;
      if (s!=s1 || s!=s2 || s!=s3 || s!=s4 || s!=s5)
         border[x] = true;
      else
         border[x] = false;
   }

   offset = xdim*(ydim-1);
   s = level[offset]>0;
   s3 = level[offset-xdim]>0;
   s4 = level[offset-xdim+1]>0;
   s5 = level[offset+1]>0;
   for (int x=1; x<xdim-1; x++)
   {
      s1 = s;
      s2 = s3;
      s = s5;
      s3 = s4;
      s4 = level[offset+x-xdim+1]>0;
      s5 = level[offset+x+1]>0;
      if (s!=s1 || s!=s2 || s!=s3 || s!=s4 || s!=s5)
         border[offset+x] = true;
      else
         border[offset+x] = false;
   }

   s = level[0]>0;
   s3 = level[1]>0;
   s4 = level[xdim+1]>0;
   s5 = level[xdim]>0;
   for (int y=1; y<ydim-1; y++)
   {
      s1 = s;
      s2 = s3;
      s = s5;
      s3 = s4;
      s4 = level[y*xdim+xdim+1]>0;
      s5 = level[y*xdim+xdim]>0;
      if (s!=s1 || s!=s2 || s!=s3 || s!=s4 || s!=s5)
         border[y*xdim] = true;
      else
         border[y*xdim] = false;
   }

   s = level[xdim-1]>0;
   s3 = level[xdim-2]>0;
   s4 = level[2*xdim-2]>0;
   s5 = level[2*xdim-1]>0;
   for (int y=1; y<ydim-1; y++)
   {
      s1 = s;
      s2 = s3;
      s = s5;
      s3 = s4;
      s4 = level[y*xdim+xdim-1+xdim-1]>0;
      s5 = level[y*xdim+xdim-1+xdim]>0;
      if (s!=s1 || s!=s2 || s!=s3 || s!=s4 || s!=s5)
         border[y*xdim+xdim-1] = true;
      else
         border[y*xdim+xdim-1] = false;
   }

   int newoffset;
   for (int y=1; y<ydim-1; y++)
   {
      offset = y*xdim;

      s4 = level[offset-xdim]>0;
      s = level[offset]>0;
      s5 = level[offset+xdim]>0;
      s6 = level[offset-xdim+1]>0;
      s7 = level[offset+1]>0;
      s8 = level[offset+xdim+1]>0;
      for (int x=1; x<xdim-1; x++)
      {
         s1 = s4;
         s2 = s;
         s3 = s5;
         s4 = s6;
         s = s7;
         s5 = s8;

         newoffset = offset + x;
         s6 = level[newoffset-xdim+1]>0;
         s7 = level[newoffset+x]>0;
         s8 = level[newoffset+xdim+1]>0;

         if (s!=s1 || s!=s2 || s!=s3 || s!=s4 || s!=s5 || s!=s6 || s!=s7 || s!=s8)
            border[newoffset] = true;
         else
            border[newoffset] = false;
      }
   }
}


