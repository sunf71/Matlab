// =============================================================================
// == NormalD.cpp
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


#include "NormalD.h"

#ifndef pi
#define pi 3.14159265
#endif

// --------------------------------------------------------------------------
// -- NormalD
// --   constructor; initializes to empty
// --------------------------------------------------------------------------
NormalD::NormalD() :
   D(0), offset(NULL), theta(NULL), Delta(NULL), total(NULL), total2(NULL), N(0),
   temp_total(NULL), temp_total2(NULL), temp_N(0), mean(NULL), uniform(false)
{
}

// --------------------------------------------------------------------------
// -- NormalD
// --   constructor; intializes to all the values given
// --------------------------------------------------------------------------
NormalD::NormalD(int newD, arr(double) newTheta, arr(double) newDelta, bool newUniform) :
   D(newD), N(0), temp_N(0), uniform(newUniform)
{
   initialize();
   sumlogDelta_div2 = 0;
   for (int d=0; d<D; d++)
   {
      offset[d] = 0;
      theta[d] = newTheta[d];
      Delta[d] = newDelta[d];
      sumlogDelta_div2 += log(Delta[d]);
      total[d] = 0;
      total2[d] = 0;
      temp_total[d] = 0;
      temp_total2[d] = 0;
   }
   sumlogDelta_div2 *= 0.5;
}
// --------------------------------------------------------------------------
// -- NormalD
// --   constructor; intializes to all the values given
// --------------------------------------------------------------------------
NormalD::NormalD(int newD, arr(double) newOffset, arr(double) newTheta, arr(double) newDelta, bool newUniform) :
   D(newD), N(0), temp_N(0), uniform(newUniform)
{
   initialize();
   sumlogDelta_div2 = 0;
   for (int d=0; d<D; d++)
   {
      offset[d] = newOffset[d];
      theta[d] = newTheta[d];
      Delta[d] = newDelta[d];
      sumlogDelta_div2 += log(Delta[d]);
      total[d] = 0;
      total2[d] = 0;
      temp_total[d] = 0;
      temp_total2[d] = 0;
   }
   sumlogDelta_div2 *= 0.5;
}
NormalD::~NormalD()
{
   cleanup();
}
void NormalD::initialize()
{
   offset = allocate_memory<double>(D);
   theta = allocate_memory<double>(D);
   Delta = allocate_memory<double>(D);
   total = allocate_memory<double>(D);
   total2 = allocate_memory<double>(D);
   temp_total = allocate_memory<double>(D);
   temp_total2 = allocate_memory<double>(D);
   mean = allocate_memory<double>(D);
}
void NormalD::cleanup()
{
   if (offset!=NULL) deallocate_memory(offset);
   if (theta!=NULL) deallocate_memory(theta);
   if (Delta!=NULL) deallocate_memory(Delta);
   if (total!=NULL) deallocate_memory(total);
   if (total2!=NULL) deallocate_memory(total2);
   if (temp_total!=NULL) deallocate_memory(temp_total);
   if (temp_total2!=NULL) deallocate_memory(temp_total2);
   if (mean!=NULL) deallocate_memory(mean);
}

// --------------------------------------------------------------------------
// -- NormalD
// --   copy constructor;
// --------------------------------------------------------------------------
NormalD::NormalD(const NormalD& that)
{
   D = that.D;
   N = that.N;
   uniform = that.uniform;
   initialize();
   for (int d=0; d<D; d++)
   {
      offset[d] = that.offset[d];
      theta[d] = that.theta[d];
      Delta[d] = that.Delta[d];
      total[d] = that.total[d];
      total2[d] = that.total2[d];
      temp_total[d] = that.temp_total[d];
      temp_total2[d] = that.temp_total2[d];
      mean[d] = that.mean[d];
   }
   sumlogDelta_div2 = that.sumlogDelta_div2;
}
// --------------------------------------------------------------------------
// -- operator=
// --   assignment operator
// --------------------------------------------------------------------------
NormalD& NormalD::operator=(const NormalD& that)
{
   if (this != &that)
   {
      cleanup();
      D = that.D;
      N = that.N;
      uniform = that.uniform;
      initialize();
      for (int d=0; d<D; d++)
      {
         offset[d] = that.offset[d];
         theta[d] = that.theta[d];
         Delta[d] = that.Delta[d];
         total[d] = that.total[d];
         total2[d] = that.total2[d];
         temp_total[d] = that.temp_total[d];
         temp_total2[d] = that.temp_total2[d];
         mean[d] = that.mean[d];
      }
      sumlogDelta_div2 = that.sumlogDelta_div2;
   }
   return *this;
}

// --------------------------------------------------------------------------
// -- copy
// --   returns a copy of this
// --------------------------------------------------------------------------
NormalD* NormalD::copy()
{
   return new NormalD(*this);
}

// --------------------------------------------------------------------------
// -- empty
// --   Empties out the NormalD. Does not update the posterior hyperparameters.
// --------------------------------------------------------------------------
void NormalD::empty()
{
   N = 0;
   temp_N = 0;

   for (int d=0; d<D; d++)
   {
      total[d] = 0;
      total2[d] = 0;
      temp_total[d] = 0;
      temp_total2[d] = 0;
   }
}

// --------------------------------------------------------------------------
// -- set_
// --   functions to setting parameters.  assumes memory already allocated
// --------------------------------------------------------------------------
void NormalD::set_offset(arr(double) newoffset)
{
   for (int d=0; d<D; d++)
      offset[d] = newoffset[d];
}
void NormalD::set_offsetsum(arr(double) offset1, arr(double) offset2)
{
   for (int d=0; d<D; d++)
      offset[d] = offset1[d] + offset2[d];
}

// --------------------------------------------------------------------------
// -- get_
// --   functions to get parameters.  returns pointers to actual data
// --------------------------------------------------------------------------
arr(double) NormalD::get_theta()
{
   if (uniform)
      for (int d=0; d<D; d++)
         theta[d] = total[d] / N;
   return theta;
}
arr(double) NormalD::get_offset()        { return offset;}
arr(double) NormalD::get_Delta()         { return Delta;}
double NormalD::get_sumlogDelta_div2()            { return sumlogDelta_div2;}
arr(double) NormalD::get_total()          { return total;}
arr(double) NormalD::get_total2()         { return total2;}

// --------------------------------------------------------------------------
// -- get_XXXX_mode
// --   returns a pointer to the most likely parameter. returns a pointer to
// -- internal memory that should *not* be deallocated
// --
// --   parameters
// --     - update : whether or not to update the posterior hyperparameters
//--    before returning the mode
// --------------------------------------------------------------------------
arr(double) NormalD::get_mean_mode()
{
   if (uniform)
      for (int d=0; d<D; d++)
         mean[d] = total[d]/N;
   else
      for (int d=0; d<D; d++)
         mean[d] = (theta[d] + offset[d] + Delta[d]*total[d]) / (1 + N*Delta[d]);
   return mean;
}
void NormalD::set_mean(arr(double) theMean)
{
   for (int d=0; d<D; d++)
      mean[d] = theMean[d];
}
void NormalD::set_meansum(arr(double) theMean1, arr(double) theMean2)
{
   for (int d=0; d<D; d++)
      mean[d] = theMean1[d] + theMean2[d];
}
arr(double) NormalD::get_mean()
{
   return mean;
}


// --------------------------------------------------------------------------
// -- add_data
// --   functions to add an observation to the NormalD. Updates the sufficient
// -- statistics stored in the data structure, but not the posterior
// -- hyperparameters
// --
// --   parameters:
// --     - data : the new observed data point of size [1 D]
// --------------------------------------------------------------------------
void NormalD::add_data(arr(double) data)
{
   // update the sufficient stats and the N
   N++;
   for (int d=0; d<D; d++)
   {
      double this_data = data[d];
      total[d] += this_data;
      total2[d] += this_data*this_data;
   }
}
void NormalD::merge_with(NormalD* &other)
{
   N += other->N;
   for (int d=0; d<D; d++)
   {
      total[d] += other->total[d];
      total2[d] += other->total2[d];
   }
}

// --------------------------------------------------------------------------
// -- rem_data
// --   functions to remove an observation from the NormalD. Updates the
// -- sufficient statistics stored in the data structure, but not the
// -- posterior hyperparameters
// --
// --   parameters:
// --     - data : the observed data point to remove of size [1 D]
// --------------------------------------------------------------------------
void NormalD::rem_data(arr(double) data)
{
   if (N<=0)
      mexErrMsgTxt("Removing a data point when there are none!\n");
   // update the sufficient stats and the N
   N--;
   for (int d=0; d<D; d++)
   {
      double this_data = data[d];
      total[d] -= this_data;
      total2[d] -= this_data*this_data;
   }
}



// --------------------------------------------------------------------------
// -- update_stats_new
// --   updates the internal posterior hyperparameters based on the stored
// -- sufficient statistics and the new data points
// --
// --   parameters:
// --     - data: the new point to consider
// --------------------------------------------------------------------------
void NormalD::update_stats_new(arr(double) data)
{
   temp_N = N + 1;
   for (int d=0; d<D; d++)
   {
      double this_data = data[d];
      temp_total[d] = total[d] + this_data;
      temp_total2[d] = total2[d] + this_data*this_data;
   }
}
void NormalD::update_stats_rem(arr(double) data)
{
   temp_N = N - 1;
   for (int d=0; d<D; d++)
   {
      double this_data = data[d];
      temp_total[d] = total[d] - this_data;
      temp_total2[d] = total2[d] - this_data*this_data;
   }
}
// --------------------------------------------------------------------------
// -- update_stats_new
// --   updates the internal posterior hyperparameters based on the stored
// -- sufficient statistics and the new data points
// --
// --   parameters:
// --     - other : another NormalD to merge with and test the likelihood
// --------------------------------------------------------------------------
void NormalD::update_stats_new(NormalD* &other)
{
   temp_N = N + other->N;
   for (int d=0; d<D; d++)
   {
      temp_total[d] = total[d] + other->total[d];
      temp_total2[d] = total2[d] + other->total2[d];
   }
}
void NormalD::update_stats_new(NormalD* &other1, NormalD* &other2)
{
   temp_N = N + other1->N + other2->N;
   for (int d=0; d<D; d++)
   {
      temp_total[d] = total[d] + other1->total[d] + other2->total[d];
      temp_total2[d] = total2[d] + other1->total2[d] + other2->total2[d];
   }
}


// --------------------------------------------------------------------------
// -- calc_loglikelihood_internal
// --   calculates log p(x | params)
// --------------------------------------------------------------------------
double NormalD::calc_logposterior_MM(arr(double) data)
{
   double logprob = -D*0.5*1.837877066409;
   for (int d=0; d<D; d++)
   {
      double temp = data[d] - mean[d];
      logprob += -0.5*(temp*temp);
   }
   return logprob;
}


// --------------------------------------------------------------------------
// -- calc_logposterior_internal
// --   calculate log p(mu|theta)p(x|mu) with either the optimal iid mu, or
// -- the current value of mu depending on fixedMean
// --------------------------------------------------------------------------
double NormalD::calc_logposterior_internal(bool useTemp, bool fixedMean)
{
   arr(double) this_total = (useTemp ? temp_total : total);
   arr(double) this_total2 = (useTemp ? temp_total2 : total2);
   double this_N = (useTemp ? temp_N : N);
   double logprob;

   if (fixedMean)
   {
      mexPrintf("NormalD shoudln't be here\n");
      // p(mu|theta)
      logprob = -D*0.5*1.837877066409 - sumlogDelta_div2;

      // p(x|mu)
      logprob += -(this_N*D)*0.5*1.837877066409;
      for (int d=0; d<D; d++)
      {
         double meand = mean[d];
         // p(mu|theta)
         double diff = meand - theta[d] - offset[d];
         logprob += -(diff*diff) / (2*Delta[d]);
         // p(x|mu)
         logprob += meand*this_total[d] - this_N*0.5*meand*meand - this_total2[d]/2;
      }
   }
   else
   {
      logprob = -D*(this_N+0.5)*0.5*1.837877066409 - sumlogDelta_div2;
      if (this_N>0)
      {
         if (uniform)
            for (int d=0; d<D; d++)
               logprob += (this_total[d]*this_total[d] / this_N - this_total2[d])/2;
         else
            for (int d=0; d<D; d++)
            {
               double totald = this_total[d];
               double totald_sq = totald*totald;
               double thetad = theta[d] + offset[d];
               logprob += (Delta[d]*totald_sq + 2*totald*thetad - this_N*thetad*thetad) / (2*this_N*Delta[d]+2) - this_total2[d]/2;
            }
      }
   }
   return logprob;
}
// --------------------------------------------------------------------------
// -- calc_logposterior
// --   calculate log p(params | x, hyperparams), for the maximal params
// -- with. Because of conjugacy, this is in the same class as
// -- log p(params | hyperparams)
// --
// --   parameters:
// --     - data : the new point to possibly add
// --------------------------------------------------------------------------
double NormalD::calc_logposterior(bool useTemp, bool fixedMean)
{
   return calc_logposterior_internal(useTemp, fixedMean);
}
// --------------------------------------------------------------------------
// -- calc_logposterior_new
// --   calculate log p(params | x, hyperparams), for the maximal params
// -- with a new data point. Because of conjugacy, this is in the same class
// -- as log p(params | hyperparams)
// --
// --   parameters:
// --     - data : the new point to possibly add
// --------------------------------------------------------------------------
double NormalD::calc_logposterior_new(arr(double) data, bool fixedMean)
{
   update_stats_new(data);
   return calc_logposterior_internal(true, fixedMean);
}
double NormalD::calc_logposterior_rem(arr(double) data, bool fixedMean)
{
   update_stats_rem(data);
   return calc_logposterior_internal(true, fixedMean);
}
// --------------------------------------------------------------------------
// -- calc_logposterior_new
// --   calculate log p(params | x, hyperparams), for the maximal params
// -- with a new merge. Because of conjugacy, this is in the same class
// -- as log p(params | hyperparams)
// --
// --   parameters:
// --     - other : another NormalD to merge with
// --------------------------------------------------------------------------
double NormalD::calc_logposterior_new(NormalD* &other, bool fixedMean)
{
   update_stats_new(other);
   return calc_logposterior_internal(true, fixedMean);
}
double NormalD::calc_logposterior_new(NormalD* &other1, NormalD* &other2, bool fixedMean)
{
   update_stats_new(other1, other2);
   return calc_logposterior_internal(true, fixedMean);
}




void NormalD::switch_priors(NormalD* &other)
{
   arr(double) temp;
   temp = offset;
   offset = other->offset;
   other->offset = temp;

   temp = theta;
   theta = other->theta;
   other->theta = temp;

   temp = Delta;
   Delta = other->Delta;
   other->Delta = temp;

   double tempDouble = sumlogDelta_div2;
   sumlogDelta_div2 = other->sumlogDelta_div2;
   other->sumlogDelta_div2 = tempDouble;

   bool tempBool = uniform;
   uniform = other->uniform;
   other->uniform = tempBool;
}

