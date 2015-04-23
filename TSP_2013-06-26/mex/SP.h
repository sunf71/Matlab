// =============================================================================
// == SP.h
// == --------------------------------------------------------------------------
// == A superpixel class to be used with temporal superpixels.
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================


#ifndef _SP_H_INCLUDED_
#define _SP_H_INCLUDED_
#include "linkedList.cpp"
#include "array.h"
#include "matrix.h"
#include "mex.h"
#include <math.h>
#include "NormalD.h"


#ifndef pi
#define pi 3.14159265
#endif

#include "helperMEX.h"
#include "gsl/gsl_sf_gamma.h"

double calc_log_label(double N, double alpha, double epsilon);
std::pair<int, int> increment_neighbor_count(std::pair<int, int> n);
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n);
std::pair<int, int> increment_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2);
std::pair<int, int> decrement_neighbor_count(std::pair<int, int> n, std::pair<int, int> n2);

class SP
{
private:
   // Likelihood term
   NormalD* pos;
   NormalD* app;

   // probabilities
   double log_likelihood;
   double log_likelihood_empty;

   // SP things
   int N;
   linkedList<int> pixels;
   linkedList<int> borders;
   linkedList< std::pair<int, int> > neighbors;

   arr(double) prev_v;

   bool is_old;

   unsigned long UID;
   friend class IMG;

public:
   // --------------------------------------------------------------------------
   // -- SP
   // --   constructor; initializes to empty... probably shouldn't use
   // --------------------------------------------------------------------------
   SP();

   // --------------------------------------------------------------------------
   // -- initialize
   // --   initializes all the variables conditioned on pos and app being set
   // --------------------------------------------------------------------------
   void initialize();
   // --------------------------------------------------------------------------
   // -- SP
   // --   constructor; initializes to empty new super pixel
   // --------------------------------------------------------------------------
   SP(NormalD &new_pos, NormalD &new_app, unsigned long new_UID, bool isOld=false, arr(double) the_prev_v=NULL);
   // --------------------------------------------------------------------------
   // -- SP
   // --   copy constructor
   // --------------------------------------------------------------------------
   SP(const SP& that);

   // --------------------------------------------------------------------------
   // -- operator=
   // --   assignment operator
   // --------------------------------------------------------------------------
   SP& operator=(const SP& that);
   // --------------------------------------------------------------------------
   // -- SP
   // --   destructor;
   // --------------------------------------------------------------------------
   ~SP();

   // --------------------------------------------------------------------------
   // -- empty
   // --   Empties out the SP
   // --
   // --   parameters:
   // --     - check_merged : checks to see if the linked list, pixels, is
   // --       empty. if it isn't, throws exception.
   // --------------------------------------------------------------------------
   void empty(bool check_merged=true);
   
   // --------------------------------------------------------------------------
   // -- isempty
   // --   returns whether or not the super pixel contains any pixels
   // --------------------------------------------------------------------------
   bool isempty();
   // --------------------------------------------------------------------------
   // -- isold
   // --   returns whether or not the superpixel is old
   // --------------------------------------------------------------------------
   bool isold();
   // --------------------------------------------------------------------------
   // -- has_no_neighbors
   // --   returns true if this SP has no neighbors
   // --------------------------------------------------------------------------
   bool has_no_neighbors();
   
   // --------------------------------------------------------------------------
   // -- checkCount
   // --   debugging function to make sure N and the linked list match in length
   // --------------------------------------------------------------------------
   void checkCount();

   // --------------------------------------------------------------------------
   // -- calculate_log_probs
   // --   Updates log_likelihood based on the current parameters
   // --------------------------------------------------------------------------
   void calculate_log_probs();

   // --------------------------------------------------------------------------
   // -- add_pixel_init
   // --   Adds a pixel to the super pixel, updates the appearance and position
   // -- parameters, and the linked lists. Does not update the likelihoods.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - is_border : indicator as to whether or not to add to border LL
   // --     - doApp (true) : indicates if the appearance should be added also
   // --   return parameters:
   // --     - pixel_ptr : a pointer to the added linked list node pixel
   // --     - border_ptr : a pointer to the added linked list node border
   // --------------------------------------------------------------------------
   void add_pixel_init(arr(double) data, int index, bool is_border,
      linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp=true);
   // --------------------------------------------------------------------------
   // -- add_pixel
   // --   Adds a pixel to the super pixel, updates the appearance and position
   // -- parameters, the linked lists, and the likelihoods. Same as
   // -- add_pixel_init except it also updates the likelihood.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - is_border : indicator as to whether or not to add to border LL
   // --     - doApp (true) : indicates if the appearance should be added also
   // --   return parameters:
   // --     - pixel_ptr : a pointer to the added linked list node pixel
   // --     - border_ptr : a pointer to the added linked list node border
   // --------------------------------------------------------------------------
   void add_pixel(arr(double) data, int index, bool is_border,
      linkedListNode<int>* &pixel_ptr, linkedListNode<int>* &border_ptr, bool doApp=true);

   // --------------------------------------------------------------------------
   // -- merge_with
   // --   Merges this SP with another one and empties the other one. Assumes
   // -- the labels have already been updated correctly. Fixes borders.
   // --
   // --   parameters:
   // --     - other : A pointer to the other SP to merge with
   // --     - label : A pointer to the label image
   // --     - border_ptr : the linked list borders to fix
   // --     - (xdim, ydim) : image dimensions
   // --------------------------------------------------------------------------
   void merge_with(SP *other, arr(int) label, arr(linkedListNode<int>*) border_ptr,
      int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- rem_pixel
   // --   Removes a pixel from the super pixel, updates the appearance and
   // -- position parameters, linked lists, and likelihoods.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - pixel_ptr : a poitner into the pixels linked list to be removed
   // --     - border_ptr : a pointer into the borders linked list to be removed
   // --     - doApp (true) : indicates if the appearance should be done also
   // --------------------------------------------------------------------------
   void rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr,
      linkedListNode<int>* &border_ptr, bool doApp=true);
   // --------------------------------------------------------------------------
   // -- rem_pixel
   // --   Removes a pixel from the super pixel, updates the appearance and
   // -- position parameters, linked lists, and likelihoods. Same as above, but
   // -- does not update any border linked lists.
   // --
   // --   parameters:
   // --     - data : [5,N] matrix containing all the data in an image
   // --     - index : in [0,N-1], index to a [5,1] data vector
   // --     - pixel_ptr : a poitner into the pixels linked list to be removed
   // --     - doApp (true) : indicates if the appearance should be done also
   // --------------------------------------------------------------------------
   void rem_pixel(arr(double) data, int index, linkedListNode<int>* &pixel_ptr, bool doApp=true);


   // --------------------------------------------------------------------------
   // -- fix_borders
   // --   Fixes the border linked list and the border_ptr image for a single
   // -- super pixel.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - border_ptr : the border_ptr image
   // --     - (xdim, ydim) : the size of the image
   // --------------------------------------------------------------------------
   void fix_borders(arr(int) label, arr(linkedListNode<int>*) border_ptr, int xdim, int ydim);



   // --------------------------------------------------------------------------
   // -- update_neighbors_label_rem
   // --   Decrements the corresponding neighbor count for neighbor_label
   // --
   // --   parameters:
   // --     - neighbor_label : the label of the neighbor to decrement
   // --------------------------------------------------------------------------
   void update_neighbors_label_rem(int neighbor_label);

   // --------------------------------------------------------------------------
   // -- update_neighbors_label_rem_check
   // --   Checks to see if the neighbor count at index should be decremented
   // -- by removing one neighbor of label neighbor_label. If so, it decrements.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index bordering the removed pixel
   // --     - (xdim,ydim) : dimensions of image
   // --     - neighbor_label : the label of the neighbor to decrement
   // --------------------------------------------------------------------------
   void update_neighbors_label_rem_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label);

   // --------------------------------------------------------------------------
   // -- update_neighbors_add_self
   // --   Updates the neighbor lists and counts by adding one particular pixel
   // -- at index. Does not update the neighboring neighbor lists.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index of the added pixel
   // --     - (xdim,ydim) : dimensions of image
   // --------------------------------------------------------------------------
   void update_neighbors_add_self(arr (int) label, int index, int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- update_neighbors_self
   // --   Updates the neighbor lists and counts by looking at all borders.
   // -- Empties previous list. The borders list must be correct!
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - (xdim,ydim) : dimensions of image
   // --------------------------------------------------------------------------
   void update_neighbors_self(arr (int) label, int xdim, int ydim);

   // --------------------------------------------------------------------------
   // -- update_neighbors_label_add
   // --   Increments the corresponding neighbor count for neighbor_label
   // --
   // --   parameters:
   // --     - neighbor_label : the label of the neighbor to increment
   // --------------------------------------------------------------------------
   void update_neighbors_label_add(int neighbor_label);
   // --------------------------------------------------------------------------
   // -- update_neighbors_label_add_check
   // --   Checks to see if the neighbor count at index should be incremented
   // -- by adding one neighbor of label neighbor_label. If so, it increments.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - label : the label image
   // --     - index : the index bordering the added pixel
   // --     - (xdim,ydim) : dimensions of image
   // --     - neighbor_label : the label of the neighbor to increment
   // --------------------------------------------------------------------------
   void update_neighbors_label_add_check(arr(int) label, int index, int xdim, int ydim, int neighbor_label);

   // --------------------------------------------------------------------------
   // -- log_likelihood_test_pointXXXX
   // --   finds the log likelihood for adding or removing a data point. All
   // -- functions take one parameter, a 5-d double vector of data. Some
   // -- functions indicated by * take an additional boolean that indicates if
   // -- the position should be checked, or the position and the appearance.
   // -- Defaults to true, which checks both.
   // -- 
   // -- XXXX can be any of the following
   // --
   // -- <empty>* : calculate the likelihood for adding
   // -- rem* : calculate the likelihood for removing
   // -- app : calculates the appearance likelihood for adding
   // -- pos : calculates the position likelihood for adding
   // -- app_rem : calculates the appearance likelihood for removing
   // -- pos_rem : calculates the position likelihood for removing
   // -- MM : appearance and position likelihood for adding w/out finding the
   // --   corresponding optimal parameters, similar to an iterative scheme
   // -- MM_pos : same as above, except only for position
   // -- 
   // --   parameters
   // --     - data : a [1 5] vector of a test point to add
   // --------------------------------------------------------------------------
   double log_likelihood_test_point(arr(double) data, bool checkApp=true);
   double log_likelihood_test_point_rem(arr(double) data, bool checkApp=true);
   double log_likelihood_test_point_app(arr(double) data);
   double log_likelihood_test_point_pos(arr(double) data);
   double log_likelihood_test_point_app_rem(arr(double) data);
   double log_likelihood_test_point_pos_rem(arr(double) data);
   double log_likelihood_test_point_MM(arr(double) data);
   double log_likelihood_test_point_MM_pos(arr(double) data);

   // --------------------------------------------------------------------------
   // -- switch_priors
   // --   switches the priors between this SP and the other SP. Neighbors
   // -- probably have to be fixed after this.
   // --
   // --   parameters:
   // --     - other : a pointer to the SP to switch priors with
   // --------------------------------------------------------------------------
   void switch_priors(SP* other);

   // --------------------------------------------------------------------------
   // -- log_likelihood_switch_app_prior
   // --   Calculates the likelihood for this SP if the appearance prior
   // -- is switched to the prior of the supplied SP
   // --
   // --   parameters:
   // --     - new_prior: a pointer to the SP to switch priors with
   // --------------------------------------------------------------------------
   double log_likelihood_switch_app_prior(SP* new_prior);
   // --------------------------------------------------------------------------
   // -- log_likelihood_switch_pos_prior
   // --   Calculates the likelihood for this SP if the position prior
   // -- is switched to the prior of the supplied SP
   // --
   // --   parameters:
   // --     - new_prior: a pointer to the SP to switch priors with
   // --------------------------------------------------------------------------
   double log_likelihood_switch_pos_prior(SP* new_prior);
   // --------------------------------------------------------------------------
   // -- log_likelihood_switch_prior
   // --   Calculates the likelihood for this SP if the prior is switched to
   // -- the prior of the supplied SP
   // --
   // --   parameters:
   // --     - new_prior: a pointer to the SP to switch priors with
   // --------------------------------------------------------------------------
   double log_likelihood_switch_prior(SP* new_prior);


   // --------------------------------------------------------------------------
   // -- log_likelihood_test_merge
   // --   Calculate the log likelihood for merging with another SP. Can be 
   // -- called with two SPs as arguments, in which case, it checks the merge
   // -- of all three SPs.  Can also use _pos to check only the position terms.
   // --
   // --   parameters
   // --     - other : another SP that we are testing for a merge
   // --------------------------------------------------------------------------
   double log_likelihood_test_merge(SP *other);
   double log_likelihood_test_merge_pos(SP *other);
   double log_likelihood_test_merge(SP *other1, SP *other2);

   // --------------------------------------------------------------------------
   // -- get_XXXX
   // --   returns a parameter.  XXXX be any of the following:
   // --
   // -- flow - the flow for this superpixel
   // -- theta_pos - the mean of the position mean
   // -- theta_app - the mean of the appearance mean
   // -- Delta_pos - the variance of the position mean
   // -- Delta_app - the variance of the appearance mean
   // -- mean_pos - the position mean
   // -- mean_app - the appearance mean
   // -- log_likelihod_pos - the log likelihood for the position params
   // -- log_likelihod_app - the log likelihood for the appearance params
   // -- prev_v - the previous flow
   // -- total_pos - the sum of all position parameters
   // -- total_app - the sum of all appearance parameters
   // -- total2_pos - the sum of all outer produces of position parameters
   // -- sumlogDelta_div2_pos - log (det(Delta)/2)
   // -- N - the number of pixels in the superpixel
   // -- log_likelihood_empty - the log likelihood if the SP was empty
   // -- log_likelihood - the log likelihood
   // -- UID - the unique ID
   // --------------------------------------------------------------------------
   arr(double) get_flow();
   arr(double) get_theta_pos();
   arr(double) get_theta_app();
   arr(double) get_Delta_pos();
   arr(double) get_Delta_app();
   arr(double) get_Sigma_pos();
   arr(double) get_Sigma_app();
   arr(double) get_mean_pos();
   arr(double) get_mean_app();
   double get_log_likelihood_app();
   double get_log_likelihood_pos();
   arr(double) get_total_app();
   arr(double) get_total_pos();
   arr(double) get_total2_pos();
   double get_sumlogDelta_div2_pos();
   int get_N();
   double get_log_likelihood_empty();
   double get_log_likelihood();
   unsigned long get_UID();
   arr(double) get_prev_v();

   // --------------------------------------------------------------------------
   // -- set_XXXX
   // --   sets the parameters indicated. XXXX be any of the following:
   // --
   // -- mean_pos - the position mean
   // -- mean_app - the appearance mean
   // -- meansum_pos - two position means to sum
   // -- flowsum - two positions to sum to the flow
   // -- flow - the flow for this superpixel (2 options)
   // --------------------------------------------------------------------------
   void set_mean_app(arr(double) new_mean);
   void set_mean_pos(arr(double) new_mean);
   void set_meansum_pos(arr(double) new_mean1, arr(double) new_mean2);
   void set_flowsum(arr(double) flow, arr(double) other);
   void set_flow(arr(double) flow);
   void set_flow(double flowx, double flowy);
};

#endif
