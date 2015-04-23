// =============================================================================
// == distancePixel.cpp
// == --------------------------------------------------------------------------
// == A class that holds a pair of values and a distance associated with the
// ==   pair.  Used for level set functions.
// == --------------------------------------------------------------------------
// == Written by Jason Chang 02-16-2008
// =============================================================================

#include "distancePixel.h"

// --------------------------------------------------------------------------
// -- pair
// --   constructor; initializes the coordinates to nothing
// --------------------------------------------------------------------------
pair::pair() :
   xpix(0), ypix(0)
{
}

// --------------------------------------------------------------------------
// -- pair
// --   constructor; initializes the coordinates to the parameters
// --------------------------------------------------------------------------   
pair::pair(const unsigned int the_xpix, const unsigned int the_ypix) :
   xpix(the_xpix), ypix(the_ypix)
{
}

// --------------------------------------------------------------------------
// -- pair
// --   copy constructor;
// --------------------------------------------------------------------------   
pair::pair(const pair &p)
{
   xpix = p.xpix;
   ypix = p.ypix;
}

// --------------------------------------------------------------------------
// -- pair
// --   assignment operator
// --------------------------------------------------------------------------
pair& pair::operator=(const pair& p)
{
   if (this != &p)
   {
      xpix = p.xpix;
      ypix = p.ypix;
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- distancePixel
// --   constructor; initializes the distancePixel to nothing
// --------------------------------------------------------------------------
distancePixel::distancePixel() :
   distance(0), coords(0,0)
{
}

// --------------------------------------------------------------------------
// -- distancePixel
// --   constructor; initializes the distancePixel to the parameters
// --------------------------------------------------------------------------
distancePixel::distancePixel(double the_distance, pair the_coords) :
   distance(the_distance), coords(the_coords)
{
}

// --------------------------------------------------------------------------
// -- distancePixel
// --   constructor; initializes the distancePixel to the parameters
// --------------------------------------------------------------------------
distancePixel::distancePixel(double the_distance, unsigned int the_xpix,
                             unsigned int the_ypix) :
   distance(the_distance), coords(the_xpix, the_ypix)
{
}

// --------------------------------------------------------------------------
// -- distancePixel
// --   assignment operator
// --------------------------------------------------------------------------
distancePixel& distancePixel::operator=(const distancePixel& p)
{
   if (this != &p)
   {
      this->distance = p.distance;
      this->coords = p.coords;
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- overloaded operators
// --   define the operators to only compare / print the distance
// --------------------------------------------------------------------------
bool distancePixel::operator==(const distancePixel dp) const
{
   return (fabs(distance) == fabs(dp.distance));
}
bool distancePixel::operator<(const distancePixel dp) const
{
   return (fabs(distance) < fabs(dp.distance));
}
bool distancePixel::operator<=(const distancePixel dp) const
{
   return (fabs(distance) <= fabs(dp.distance));
}
bool distancePixel::operator>(const distancePixel dp) const
{
   return (fabs(distance) > fabs(dp.distance));
}
bool distancePixel::operator>=(const distancePixel dp) const
{
   return (fabs(distance) >= fabs(dp.distance));
}
std::ostream& operator<< (std::ostream& os, distancePixel dp)
{
   os << dp.distance << " (" << dp.coords.xpix<< "," << dp.coords.ypix<<")";
   return os;
}

