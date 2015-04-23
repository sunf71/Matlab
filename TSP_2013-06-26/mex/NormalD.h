// =============================================================================
// == NormalD.h
// == --------------------------------------------------------------------------
// == A class of D-dimensional Normal random variables with unknown mean and
// == fixed variance. Assumes a Normal prior on the mean.
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================


#ifndef _NormalD_H_INCLUDED_
#define _NormalD_H_INCLUDED_

#include "mex.h"
#include <math.h>
#include "array.h"
#include "linear_algebra.h"

#include "helperMEX.h"
#include "debugMEX.h"

class NormalD
{
protected:
   int D;
   // prior hyperparameters
   // x ~ N(mean, Sigma)
   // mean ~ N(theta, Delta);
   arr(double) offset;
   arr(double) theta;
   arr(double) Delta;

   double sumlogDelta_div2; // 0.5*sum(log(Delta(:)))

   // sufficient statistics of the observed data
   arr(double) total;
   arr(double) total2;
   int N;

   // temporary sufficient statistics
   arr(double) temp_total;
   arr(double) temp_total2;
   int temp_N;

   // instantiated gaussian parameters
   arr(double) mean;

   bool uniform;

public:
   // --------------------------------------------------------------------------
   // -- NormalD
   // --   constructor; initializes to empty
   // --------------------------------------------------------------------------
   NormalD();
   // --------------------------------------------------------------------------
   // -- NormalD
   // --   constructor; intializes to all the values given
   // --------------------------------------------------------------------------
   NormalD(int newD, arr(double) newTheta, arr(double) newDelta, bool newUniform=false);
   // --------------------------------------------------------------------------
   // -- NormalD
   // --   constructor; intializes to all the values given
   // --------------------------------------------------------------------------
   NormalD(int newD, arr(double) newOffset, arr(double) newTheta, arr(double) newDelta, bool newUniform=false);
   ~NormalD();
   void initialize();
   void cleanup();

   // --------------------------------------------------------------------------
   // -- NormalD
   // --   copy constructor;
   // --------------------------------------------------------------------------
   NormalD(const NormalD& that);
   // --------------------------------------------------------------------------
   // -- operator=
   // --   assignment operator
   // --------------------------------------------------------------------------
   NormalD& operator=(const NormalD& that);
   // --------------------------------------------------------------------------
   // -- copy
   // --   returns a copy of this
   // --------------------------------------------------------------------------
   NormalD* copy();

   // --------------------------------------------------------------------------
   // -- empty
   // --   Empties out the NormalD. Does not update the posterior hyperparameters.
   // --------------------------------------------------------------------------
   void empty();

   // --------------------------------------------------------------------------
   // -- set_
   // --   functions to setting parameters.  assumes memory already allocated
   // --------------------------------------------------------------------------
   void set_offset(arr(double) newoffset);
   void set_offsetsum(arr(double) offset1, arr(double) offset2);

   // --------------------------------------------------------------------------
   // -- get_
   // --   functions to get parameters.  returns pointers to actual data
   // --------------------------------------------------------------------------
   arr(double) get_theta();
   arr(double) get_offset();
   arr(double) get_Delta();
   arr(double) get_total();
   arr(double) get_total2();
   double get_sumlogDelta_div2();

   // --------------------------------------------------------------------------
   // -- get_XXXX_mode
   // --   returns a pointer to the most likely parameter. returns a pointer to
   // -- internal memory that should *not* be deallocated
   // --
   // --   parameters
   // --     - update : whether or not to update the posterior hyperparameters
   //--    before returning the mode
   // --------------------------------------------------------------------------
   arr(double) get_mean_mode();
   void set_mean(arr(double) theMean);
   void set_meansum(arr(double) theMean1, arr(double) theMean2);
   arr(double) get_mean();

   // --------------------------------------------------------------------------
   // -- add_data
   // --   functions to add an observation to the NormalD. Updates the sufficient
   // -- statistics stored in the data structure, but not the posterior
   // -- hyperparameters
   // --
   // --   parameters:
   // --     - data : the new observed data point of size [1 D]
   // --------------------------------------------------------------------------
   void add_data(arr(double) data);

   void merge_with(NormalD* &other);
   // --------------------------------------------------------------------------
   // -- rem_data
   // --   functions to remove an observation from the NormalD. Updates the
   // -- sufficient statistics stored in the data structure, but not the
   // -- posterior hyperparameters
   // --
   // --   parameters:
   // --     - data : the observed data point to remove of size [1 D]
   // --------------------------------------------------------------------------
   void rem_data(arr(double) data);



   double calc_logposterior_MM(arr(double) data);

   // --------------------------------------------------------------------------
   // -- calc_logposterior
   // --   calculate log p(params | x, hyperparams), for the maximal params
   // -- with. Because of conjugacy, this is in the same class as
   // -- log p(params | hyperparams)
   // --
   // --   parameters:
   // --     - data : the new point to possibly add
   // --------------------------------------------------------------------------
   double calc_logposterior(bool useTemp=false, bool fixedMean=false);
   // --------------------------------------------------------------------------
   // -- calc_logposterior_new
   // --   calculate log p(params | x, hyperparams), for the maximal params
   // -- with a new data point. Because of conjugacy, this is in the same class
   // -- as log p(params | hyperparams)
   // --
   // --   parameters:
   // --     - data : the new point to possibly add
   // --------------------------------------------------------------------------
   double calc_logposterior_new(arr(double) data, bool fixedMean=false);
   double calc_logposterior_rem(arr(double) data, bool fixedMean=false);
   // --------------------------------------------------------------------------
   // -- calc_logposterior_new
   // --   calculate log p(params | x, hyperparams), for the maximal params
   // -- with a new merge. Because of conjugacy, this is in the same class
   // -- as log p(params | hyperparams)
   // --
   // --   parameters:
   // --     - other : another NormalD to merge with
   // --------------------------------------------------------------------------
   double calc_logposterior_new(NormalD* &other, bool fixedMean=false);
   double calc_logposterior_new(NormalD* &other1, NormalD* &other2, bool fixedMean=false);


   void switch_priors(NormalD* &other);


private:
   // --------------------------------------------------------------------------
   // -- update_posteriors_new
   // --   updates the internal posterior hyperparameters based on the stored
   // -- sufficient statistics and the new data points
   // --
   // --   parameters:
   // --     - data: the new point to consider
   // --------------------------------------------------------------------------
   void update_stats_new(arr(double) data);
   void update_stats_rem(arr(double) data);
   // --------------------------------------------------------------------------
   // -- update_posteriors_new
   // --   updates the internal posterior hyperparameters based on the stored
   // -- sufficient statistics and the new data points
   // --
   // --   parameters:
   // --     - other : another NormalD to merge with and test the likelihood
   // --------------------------------------------------------------------------
   void update_stats_new(NormalD* &other);
   void update_stats_new(NormalD* &other1, NormalD* &other2);

   // --------------------------------------------------------------------------
   // -- calc_logposterior_internal
   // --   calculate log p(mu|theta)p(x|mu) with either the optimal iid mu, or
   // -- the current value of mu depending on fixedMean
   // --------------------------------------------------------------------------
   double calc_logposterior_internal(bool useTemp=false, bool fixedMean=false);
};


#endif
