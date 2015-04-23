// =============================================================================
// == IMG.h
// == --------------------------------------------------------------------------
// == An image class to be used with temporal superpixels.
// ==
// == All work using this code should cite:
// == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
// ==    Temporal Superpixels. CVPR 2013.
// == --------------------------------------------------------------------------
// == Written by Jason Chang and Donglai Wei 06-20-2013
// =============================================================================

#ifndef _IMG_H_INCLUDED_
#define _IMG_H_INCLUDED_
#include <iostream>
#include <fstream>
#include "utils.h"
#include "SP.h"
#include "NormalD.h"
#include "topology.h"
#include "array.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "helperMEX.h"

using namespace std;

class IMG
{
private:
   int debug;
   int xdim;
   int ydim;
   int N;

   arr(double) data;
   arr(int) label;
   int K;
   //NIW new_pos;
   NormalD new_pos;
   NormalD new_app;
   arr(SP*) SP_arr;
   SP* SP_new; // a pointer to an exemplar "new" SP
   arr(bool) SP_old; // indicates the SPs that are from a prev frame
   arr(bool) SP_changed; // indicates which SPs changed in the previous iteration

   arr(linkedListNode<int>*) border_ptr;
   arr(linkedListNode<int>*) pixel_ptr;

   arr(double) T4Table;

   double log_alpha; // size parameter
   double log_beta;
   double area;
   double area_var;
   double log_area_var;
   unsigned long max_UID;

   double lambda; // flow parameters
   double lambda_sigma;

   arr(int) prev_label;
   arr(double) prev_app_mean;
   arr(double) prev_pos_mean;
   int prev_K;
   arr(double) prev_covariance;
   arr(double) prev_precision;
   linkedList< std::pair<int, int> > prev_neighbors;
   linkedList< double > prev_neighbors_scale;
   bool alive_dead_changed;

   double dummy_log_prob;
   arr(bool) boundary_mask;

   // flow stuff
   gsl_vector *vgsl;
   gsl_vector *ugsl;
   arr(double) Sxy;
   arr(double) Syy;
   arr(double) SxySyy;
   arr(double) obs_u;
   arr(double) obs_v;
   arr(double) temp_Syy;
   arr(double) temp_d;
   arr(bool) temp_still_alive;
   int num_alive;
   int num_dead;

   arr(int) all2alive;
   arr(int) alive2all;
   arr(int) all2dead;
   arr(int) dead2all;

public:
   IMG();
   virtual ~IMG();

   void ReadIMG();
   void mxReadIMG(const mxArray *mstruct);
   void mxWriteIMG(mxArray* plhs[], const mxArray* oldstruct);


   // --------------------------------------------------------------------------
   // -- move_local
   // --   finds the optimal local joint move of labels and parameters. Chooses
   // -- a super pixel at random and loops through and updates its borders.
   // --------------------------------------------------------------------------
   bool move_local_IMG();


   bool move_switch_IMG();


   // --------------------------------------------------------------------------
   // -- move_merge
   // --   finds the optimal merge between all pairs of super pixels
   // --------------------------------------------------------------------------
   void move_merge_IMG();
   void move_split_IMG();
   void move_IMG();

   void Flow_QP2();
   void find_SxySyy(arr(int) all2alive, arr(int) alive2all, arr(int) all2dead, arr(int) dead2all);
private:
   //
   // -- 1. I/O
   void U_initialize();
   void U_relabel_SP(bool final=false);


   double U_calc_energy();
   double U_calc_model_order(int N, bool is_old);

   // -- 3. border bookkeeping
   bool U_check_border_pix(int index);
   bool U_check_border_pix(int index, int cur_label);
   void U_find_border_SP(int k,arr(bool) neighbors,linkedList<int> &neighborsLL);
   // --------------------------------------------------------------------------
   // -- U_update_border_changed
   // --   Updates all border linked lists and pointers for the neighbors of a
   // -- changed pixel at index.
   // --
   // --   parameters:
   // --     - index : the index of the recently changed pixel
   // --------------------------------------------------------------------------
   void U_update_border_changed(int index);
   // --------------------------------------------------------------------------
   // -- U_update_border_changed_pixel
   // --   Updates all border linked lists and pointers at index
   // --
   // --   parameters:
   // --     - index : the index of the recently changed pixel
   // --------------------------------------------------------------------------
   void U_update_border_changed_pixel(int index);

   // --------------------------------------------------------------------------
   // -- U_update_neighbor_list
   // --   Updates super pixel neighbors used in merge
   // --
   // --   parameters:
   // --     - neighbors : a boolean array indicating which labels are added
   // --     - neighborsLL : a linked list of the neighbor indices
   // --     - index : the index of the point to consider adding
   // --------------------------------------------------------------------------
   void U_update_neighbor_list(arr(bool) neighbors, linkedList<int> &neighborsLL, int index);


   void U_fix_neighbors_self(int k);
   void U_fix_neighbors_neighbors(int k);
   void U_fix_neighbors_neighbors(int k, int kignore);
   void U_print_neighbors(int k);

   // --------------------------------------------------------------------------
   // -- U_update_neighbors_rem
   // --   Updates all neighbor lists when a pixel is "removed". A subsequent
   // -- call to update_neighbors_add should be completed right after this one.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - old_label : the label of the pixel before it was changed
   // --     - index : the index bordering the removed pixel
   // --------------------------------------------------------------------------
   void U_update_neighbors_rem(int old_label, int index);

   // --------------------------------------------------------------------------
   // -- U_update_neighbors_add
   // --   Updates all neighbor lists when a pixel is "added". A previous
   // -- call to update_neighbors_rem should be completed right before this one.
   // -- The neighboring label should be changed before calling this function.
   // --
   // --   parameters:
   // --     - index : the index bordering the removed pixel
   // --------------------------------------------------------------------------
   void U_update_neighbors_add(int index);

   // --------------------------------------------------------------------------
   // -- U_update_neighbors_merge
   // --   Updates the neighbor lists for merging two super pixels.
   // -- All labels and SP_arr things should be updated *before* calling this.
   // --
   // --   parameters:
   // --     - index : the index bordering the removed pixel
   // --------------------------------------------------------------------------
   void U_update_neighbors_merge(int new_label, int old_label);

   // --------------------------------------------------------------------------
   // -- U_update_neighbors_split
   // --   Updates the neighbor lists for merging two super pixels.
   // -- All labels and SP_arr things should be updated *before* calling this.
   // --
   // --   parameters:
   // --     - index : the index bordering the removed pixel
   // --------------------------------------------------------------------------
   void U_update_neighbors_split(int label1, int label2);

   // -- 4. Modularized Move
   // -- 4.0 Change of Energy

   // --------------------------------------------------------------------------
   // -- link_cost
   // --   calculates the energy of linking SP_arr[ko] (old) to SP_arr[kn] (new)
   // --
   // --   parameters:
   // --     - ko : the old k
   // --     - kn : the new k
   // --------------------------------------------------------------------------
   double link_cost(int ko, int kn);

   double move_local_calc_delta_MM(int index, int new_k, double& max_prob, int& max_k);

   // --------------------------------------------------------------------------
   // -- move_local_calc_neigbor
   // --   calculates the probability of assigning the pixel at index to the
   // -- cluster of nindex (neighbor index). Updates the max_prob and max_k if
   // -- the new_app value is greater than the old one
   // --
   // --   parameters:
   // --     - index : the new point to add
   // --     - new_k : the neighbor to add this point to
   // --     - max_prob : the maximum probability of all neighbors
   // --     - max_k : the super pixel index of the maximum probability
   // --------------------------------------------------------------------------
   double move_local_calc_delta(int index, int new_k, bool add, double& max_prob, int& max_k);

   double move_switch_calc_delta(SP* &oldSP, SP* &newSP);


   // --------------------------------------------------------------------------
   // -- move_merge_calc_delta
   // --   calculates the probability of assigning the pixel at index to the
   // -- cluster of nindex (neighbor index).
   // --
   // --   parameters:
   // --     - index : the new point to add
   // --     - nindex : the neighbor to add this point to
   // --------------------------------------------------------------------------
   double move_merge_calc_delta(int k, int merge_k);


   // --------------------------------------------------------------------------
   // -- move_split_calc_delta
   // --   calculates the change in energy for
   // -- (k1 U new_k1) && (k2 U new_k2) - (k1 U new_k1 U new_k) && (k2)
   // --
   // --   parameters:
   // --     - SP1 : the SP that originates the splitting
   // --     - SP2 : the SP to split to
   // --     - new_SP1 : temporary SP that contains pixels that will go in k1
   // --     - new_SP2 : temporary SP that contains pixels that will go in k2
   // --     - SP1_old : indicates if SP1 is an old SP
   // --     - SP2_old : indicates if SP2 is an old SP
   // --------------------------------------------------------------------------
   double move_split_calc_delta(SP* SP1, SP* SP2, SP* new_SP1, SP* new_SP2, bool SP1_old, bool SP2_old);

   // -- 4.1 move local
   // need to be done
   int move_local_SP_region(int k,linkedList<int> &check_labels);
   int move_local_SP(int k);
   // need to be done
   int move_local_pix(int index);
   //-- 4.2 move merge
   void move_merge_SP_propose(int k,arr(bool) neighbors,double &max_E,int &max_k);
   void move_merge_SP_propose_region(int k,arr(bool) neighbors, linkedList<int> &check_labels,double &max_E,int &max_k);
   void move_merge_SP(int k,arr(bool) neighbors);

   // -- 4.3 move split
   void move_split_SP(int k);
   bool U_Kmeans_plusplus(int k, int *bbox, int num_SP, int numiter, int index2=-1);
   void U_connect_newSP(int *bbox, int num_SP);
   double U_dist(int index1 , double* center);
   double U_dist(int index1 , int index2);

   void move_split_SP_propose(int index, int num_SP, int option, double &max_E,int &ksplit,int* new_ks);



};

#endif
