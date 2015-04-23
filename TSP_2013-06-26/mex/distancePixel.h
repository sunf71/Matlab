// =============================================================================
// == distancePixel.h
// == --------------------------------------------------------------------------
// == A class that holds a pair of values and a distance associated with the
// ==   pair.  Used for level set functions.
// == --------------------------------------------------------------------------
// == Written by Jason Chang 02-16-2008
// =============================================================================

#ifndef distancePixel_H
#define distancePixel_H

#include <iostream>
#include <math.h>

class pair
{
private:
   friend class distancePixel;
   friend class level_set;
   friend class m_level_set;

public:
   unsigned int xpix;
   unsigned int ypix;
   // --------------------------------------------------------------------------
   // -- pair
   // --   constructor; initializes the coordinates to nothing
   // --------------------------------------------------------------------------
   pair();

   // --------------------------------------------------------------------------
   // -- pair
   // --   constructor; initializes the coordinates to the parameters
   // --------------------------------------------------------------------------   
   pair(const unsigned int the_xpix, const unsigned int the_ypix);
   
   // --------------------------------------------------------------------------
   // -- pair
   // --   copy constructor;
   // --------------------------------------------------------------------------   
   pair(const pair &p);
   
   // --------------------------------------------------------------------------
   // -- pair
   // --   assignment operator
   // --------------------------------------------------------------------------
   pair& operator=(const pair& p);
};

class distancePixel
{
private:
   friend class level_set;

public:
   double distance;
   pair coords;

   // --------------------------------------------------------------------------
   // -- distancePixel
   // --   constructor; initializes the distancePixel to nothing
   // --------------------------------------------------------------------------
   distancePixel();

   // --------------------------------------------------------------------------
   // -- distancePixel
   // --   constructor; initializes the distancePixel to the parameters
   // --------------------------------------------------------------------------
   distancePixel(double the_distance, pair the_coords);

   // --------------------------------------------------------------------------
   // -- distancePixel
   // --   constructor; initializes the distancePixel to the parameters
   // --------------------------------------------------------------------------
   distancePixel(double the_distance, unsigned int the_xpix,
                 unsigned int the_ypix);

   // --------------------------------------------------------------------------
   // -- distancePixel
   // --   assignment operator
   // --------------------------------------------------------------------------
   distancePixel& operator=(const distancePixel& p);

   // --------------------------------------------------------------------------
   // -- overloaded operators
   // --   define the operators to only compare / print the distance
   // --------------------------------------------------------------------------
   bool operator==(const distancePixel dp) const;
   bool operator<(const distancePixel dp) const;
   bool operator<=(const distancePixel dp) const;
   bool operator>(const distancePixel dp) const;
   bool operator>=(const distancePixel dp) const;
   friend std::ostream& operator<< (std::ostream& os, distancePixel dp);
};

#endif /* distancePixel_H */

