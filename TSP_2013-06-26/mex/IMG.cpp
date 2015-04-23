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


#include "IMG.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_cdf.h"

//////////////////////////////////////////////////////////////////////
// 1. Construction/Destruction
//////////////////////////////////////////////////////////////////////
IMG::IMG() : K(0)
{
	label = NULL;
	SP_arr = NULL;
	SP_new = NULL;
	SP_old = NULL;
	SP_changed = NULL;
	border_ptr = NULL;
	pixel_ptr = NULL;
   prev_label = NULL;
}

IMG::~IMG()
{
   if (SP_arr)
   {
      for (int k=0; k<N; k++)
         if (SP_arr[k]!=NULL)
            delete SP_arr[k];
      deallocate_memory(SP_arr);
   }
   if (SP_new)
      delete SP_new;
   if (SP_old)
      deallocate_memory(SP_old);
   if (SP_changed)
      deallocate_memory(SP_changed);
   if (border_ptr)
      deallocate_memory(border_ptr);
   if (pixel_ptr)
      deallocate_memory(pixel_ptr);
   if (label)
      deallocate_memory(label);

   if (vgsl != NULL)
      gsl_vector_free(vgsl);
   if (ugsl != NULL)
      gsl_vector_free(ugsl);
   if (Sxy != NULL)
      deallocate_memory(Sxy);
   if (Syy != NULL)
      deallocate_memory(Syy);
   if (SxySyy != NULL)
      deallocate_memory(SxySyy);
   if (obs_u != NULL)
      deallocate_memory(obs_u);
   if (obs_v != NULL)
      deallocate_memory(obs_v);
   if (temp_Syy != NULL)
      deallocate_memory(temp_Syy);
   if (temp_d != NULL)
      deallocate_memory(temp_d);
   if (temp_still_alive != NULL)
      deallocate_memory(temp_still_alive);

   if (all2alive != NULL)
      deallocate_memory(all2alive);
   if (alive2all != NULL)
      deallocate_memory(alive2all);
   if (all2dead != NULL)
      deallocate_memory(all2dead);
   if (dead2all != NULL)
      deallocate_memory(dead2all);
}
//////////////////////////////////////////////////////////////////////
// 2. I/O
//////////////////////////////////////////////////////////////////////
void IMG::mxReadIMG(const mxArray *mstruct)
{
   int nfields = mxGetNumberOfFields(mstruct);
   int nIMG = mxGetNumberOfElements(mstruct);
   arr(double) tmp_arr;

   /*if (nfields!=11 && nfields!=12 && nfields!=19){mexErrMsgTxt("Input IMG structure must have only 10 or 11 fields\n");}
   if (nIMG != 1){mexErrMsgTxt("Input IMG structure must have only 1 element\n");}*/

   // read in the variables from the structure
   xdim = (int)mxReadScalar(mxGetField(mstruct, 0, "xdim"));
   ydim = (int)mxReadScalar(mxGetField(mstruct, 0, "ydim"));
   N = xdim*ydim;
   K = (int)mxReadScalar(mxGetField(mstruct, 0, "K"));
   data = getArrayInput<double>(mxGetField(mstruct,0,"data"));
   T4Table = getArrayInput<double>(mxGetField(mstruct, 0, "T4Table"));

   log_alpha = (double)mxReadScalar(mxGetField(mstruct, 0, "log_alpha"));
   log_beta = (double)mxReadScalar(mxGetField(mstruct, 0, "log_beta"));
   area = (double)mxReadScalar(mxGetField(mstruct, 0, "area"));
   area_var = (double)mxReadScalar(mxGetField(mstruct, 0, "area_var"));
   log_area_var = log(area_var);
   max_UID = getInput<unsigned long>(mxGetField(mstruct, 0, "max_UID"));

   // create temporary space for variables that we will change
   checkInput(mxGetField(mstruct,0,"label"), MATRIX, "int32");
   int* temp_label = getArrayInput<int>(mxGetField(mstruct, 0, "label"));
   label = allocate_memory<int>(N);
   memcpy(label, temp_label, sizeof(int)*N);

   // create the image pointer stuff
   pixel_ptr = allocate_memory<linkedListNode<int>*>(N);
   border_ptr = allocate_memory<linkedListNode<int>*>(N);

   // read in the hyper parameters for the new super pixel
   mxArray* hyper = mxGetField(mstruct, 0, "hyper");

   new_pos = NormalD(2,
      getArrayInput<double>(mxGetField(hyper,0,"p_theta")),
      getArrayInput<double>(mxGetField(hyper,0,"p_Delta")),
      true);
   new_app = NormalD(3,
      getArrayInput<double>(mxGetField(hyper,0,"a_theta")),
      getArrayInput<double>(mxGetField(hyper,0,"a_Delta")),
      true);
   SP_new = new SP(new_pos, new_app, 0, false);

   // read in the super pixel array
   mxArray* tmp_SP = mxGetField(mstruct,0,"SP");
   int num_old = mxGetNumberOfElements(tmp_SP);

   SP_arr = allocate_memory<SP*>(N);
   SP_old = allocate_memory<bool>(N, false);

   for (int i=0; i<num_old; i++)
   {
      mxArray* tmp = mxGetField(tmp_SP,i,"old");
      if (tmp==NULL || mxIsEmpty(tmp) || ! (bool)mxGetScalar(tmp))
      {
         tmp = mxGetField(tmp_SP,i,"UID");
         if (tmp==NULL || mxIsEmpty(tmp))
            SP_arr[i] = new SP(new_pos, new_app, max_UID++, false);
         else
            SP_arr[i] = new SP(new_pos, new_app, getInput<unsigned long>(tmp), false);
         SP_old[i] = false;
      }
      else
      {
         NormalD tmp_pos = NormalD(2,
            getArrayInput<double>(mxGetField(tmp_SP,i,"v")),
            getArrayInput<double>(mxGetField(tmp_SP,i,"p_theta")),
            getArrayInput<double>(mxGetField(tmp_SP,i,"p_Delta")),
            false);
         NormalD tmp_app = NormalD(3,
            getArrayInput<double>(mxGetField(tmp_SP,i,"a_theta")),
            getArrayInput<double>(mxGetField(tmp_SP,i,"a_Delta")),
            false);
         unsigned long UID = getInput<unsigned long>(mxGetField(tmp_SP,i,"UID"));
         arr(double) prev_v = getArrayInput<double>(mxGetField(tmp_SP,i,"prev_v"));
         SP_arr[i] = new SP(tmp_pos, tmp_app, UID, true, prev_v);
         SP_arr[i]->set_mean_pos(getArrayInput<double>(mxGetField(tmp_SP,i,"p_mu")));
         SP_arr[i]->set_mean_app(getArrayInput<double>(mxGetField(tmp_SP,i,"a_mu")));
         SP_old[i] = true;
      }
      /*tmp = mxGetField(tmp_SP,i,"old");
      if (tmp!=NULL && !mxIsEmpty(tmp))
      {
         mexPrintf("%d %f\n", (*(bool*)(tmp)), mxGetScalar(tmp));
      }*/
   }
   for (int i=num_old; i<K; i++)
      SP_arr[i] = new SP(new_pos, new_app, max_UID++, false);
   for (int i=K; i<N; i++)
      SP_arr[i] = NULL;

   SP_changed = allocate_memory<bool>(N);
   mxArray* tmp = mxGetField(mstruct, 0, "SP_changed");
   if (tmp==NULL) // doesn't exist... create it
   {
      for (int i=0; i<N; i++)
         SP_changed[i] = true;
   }
   else
   {
      arr(bool) tmp_SP_changed = getArrayInput<bool>(mxGetField(mstruct,0,"SP_changed"));
      for (int i=0; i<N; i++)
         SP_changed[i] = tmp_SP_changed[i];
   }

   boundary_mask = getArrayInput<bool>(mxGetField(mstruct,0,"boundary_mask"));
   dummy_log_prob = getInput<double>(mxGetField(mstruct,0,"dummy_log_prob"));

   // read in the previous labels
   if (mxGetFieldNumber(mstruct, "prev_label")>0)
   {
      prev_label = getArrayInput<int>(mxGetField(mstruct,0,"prev_label"));
      prev_app_mean = getArrayInput<double>(mxGetField(mstruct,0,"prev_app_mean"));
      prev_pos_mean = getArrayInput<double>(mxGetField(mstruct,0,"prev_pos_mean"));
      prev_K = getInput<double>(mxGetField(mstruct,0,"prev_K"));
      //lambda = getInput<double>(mxGetField(mstruct,0,"lambda"));
      //lambda_sigma = getInput<double>(mxGetField(mstruct,0,"lambda_sigma"));
      prev_covariance = getArrayInput<double>(mxGetField(mstruct,0,"prev_covariance"));
      prev_precision = getArrayInput<double>(mxGetField(mstruct,0,"prev_precision"));
      alive_dead_changed = getInput<bool>(mxGetField(mstruct,0,"alive_dead_changed"));
      if (alive_dead_changed)
      {
         num_alive = 0;
         Sxy = NULL;
         Syy = NULL;
         SxySyy = NULL;
         obs_u = NULL;
         obs_v = NULL;
      }
      else
      {
         arr(double) tempSxy = getArrayInput<double>(mxGetField(mstruct,0,"Sxy"));
         arr(double) tempSyy = getArrayInput<double>(mxGetField(mstruct,0,"Syy"));
         arr(double) tempSxySyy = getArrayInput<double>(mxGetField(mstruct,0,"SxySyy"));
         num_alive = mxGetM(mxGetField(mstruct,0,"SxySyy"));
         Sxy = allocate_memory<double>(num_alive*prev_K);
         Syy = allocate_memory<double>(num_alive*num_alive);
         SxySyy = allocate_memory<double>(num_alive*prev_K);
         memcpy(Sxy, tempSxy, sizeof(double)*num_alive*prev_K);
         memcpy(Syy, tempSyy, sizeof(double)*num_alive*num_alive);
         memcpy(SxySyy, tempSxySyy, sizeof(double)*num_alive*prev_K);

         obs_u = allocate_memory<double>(num_alive);
         obs_v = allocate_memory<double>(num_alive);
      }

      ugsl = gsl_vector_alloc(prev_K);
      vgsl = gsl_vector_alloc(prev_K);
      temp_Syy = allocate_memory<double>(prev_K*prev_K);
      temp_d = allocate_memory<double>(prev_K);
      temp_still_alive = allocate_memory<bool>(prev_K);

      all2alive = allocate_memory<int>(prev_K);
      alive2all = allocate_memory<int>(prev_K);
      all2dead = allocate_memory<int>(prev_K);
      dead2all = allocate_memory<int>(prev_K);
   }
   else // first frame
   {
      prev_label = NULL;
      prev_app_mean = NULL;
      prev_pos_mean = NULL;
      prev_K = 0;
      //lambda = 0;
      //lambda_sigma = 0;
      alive_dead_changed = true;
      num_alive = 0;
      Sxy = NULL;
      Syy = NULL;
      SxySyy = NULL;
      obs_u = NULL;
      obs_v = NULL;
      vgsl = NULL;
      ugsl = NULL;
      temp_Syy = NULL;
      temp_d = NULL;
      temp_still_alive = NULL;

      all2alive = NULL;
      alive2all = NULL;
      all2dead = NULL;
      dead2all = NULL;
   }
   U_initialize();
   //mexPrintf("%f\n", U_calc_energy());
}

void IMG::mxWriteIMG(mxArray* plhs[], const mxArray* oldstruct)
{
   plhs[0] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
   arr(double) temp_K = getArrayInput<double>(plhs[0]);
   temp_K[0] = K;

   plhs[1] = mxCreateNumericMatrix(xdim, ydim, mxINT32_CLASS, mxREAL);
   arr(int) temp_label = getArrayInput<int>(plhs[1]);
   memcpy(temp_label, label, sizeof(int)*N);

   const char* SP_names[11] = {"p_Delta", "p_theta", "p_mu", "a_Delta", "a_theta", "a_mu", "N", "UID", "v", "old", "prev_v"};
   plhs[2] = mxCreateStructMatrix(K,1,11,SP_names);
   mxArray* tmp_SP = plhs[2];
   for (int i =0; i<K; i++)
   {
      if (SP_arr[i]!=NULL)
      {
         // Pixel stats
         mxWriteField(tmp_SP,i,"p_theta",mxWriteVector2<double>(1,2,    SP_arr[i]->get_theta_pos()));
         mxWriteField(tmp_SP,i,"v",      mxWriteVector2<double>(1,2,    SP_arr[i]->get_flow()));
         mxWriteField(tmp_SP,i,"p_Delta",mxWriteVector2<double>(1,2,    SP_arr[i]->get_Delta_pos()));
         mxWriteField(tmp_SP,i,"p_mu",   mxWriteVector2<double>(1,2,    SP_arr[i]->get_mean_pos()));
         mxWriteField(tmp_SP,i,"a_Delta",mxWriteVector2<double>(1,3,    SP_arr[i]->get_Delta_app()));
         mxWriteField(tmp_SP,i,"a_theta",mxWriteVector2<double>(1,3,    SP_arr[i]->get_theta_app()));
         mxWriteField(tmp_SP,i,"a_mu",   mxWriteVector2<double>(1,3,    SP_arr[i]->get_mean_app()));
         mxWriteField(tmp_SP,i,"N",      mxWriteScalar(        (double)SP_arr[i]->get_N()));
         mxWriteField(tmp_SP,i,"UID",    mxWriteScalar(                SP_arr[i]->get_UID()));
         mxWriteField(tmp_SP,i,"old",    mxWriteScalar(                SP_old[i]));
         mxWriteField(tmp_SP,i,"prev_v", mxWriteVector2<double>(1,2,    SP_arr[i]->get_prev_v()));
      }
      else if (SP_old[i])
      {
         mexPrintf("AHHHHHHHHHH!\n");
      }
   }

   plhs[3] = mxCreateLogicalMatrix(1, N);
   arr(bool) temp_SP_changed = getArrayInput<bool>(plhs[3]);
   for (int i=0; i<N; i++)
      temp_SP_changed[i] = SP_changed[i] && SP_arr[i]!=NULL && !SP_arr[i]->isempty();

   plhs[4] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
   arr(unsigned long) temp_max_UID = getArrayInput<unsigned long>(plhs[4]);
   temp_max_UID[0] = max_UID;

   plhs[5] = mxCreateLogicalMatrix(1, 1);
   arr(bool) temp_alive_dead_changed = getArrayInput<bool>(plhs[5]);
   temp_alive_dead_changed[0] = alive_dead_changed;

   plhs[6] = mxCreateNumericMatrix(num_alive, prev_K, mxDOUBLE_CLASS, mxREAL);
   arr(double) tempSxy = getArrayInput<double>(plhs[6]);
   memcpy(tempSxy, Sxy, sizeof(double)*num_alive*prev_K);

   plhs[7] = mxCreateNumericMatrix(num_alive, num_alive, mxDOUBLE_CLASS, mxREAL);
   arr(double) tempSyy = getArrayInput<double>(plhs[7]);
   memcpy(tempSyy, Syy, sizeof(double)*num_alive*num_alive);

   plhs[8] = mxCreateNumericMatrix(num_alive, prev_K, mxDOUBLE_CLASS, mxREAL);
   arr(double) tempSxySyy = getArrayInput<double>(plhs[8]);
   memcpy(tempSxySyy, SxySyy, sizeof(double)*num_alive*prev_K);

   //mexPrintf("%f\n", U_calc_energy());
   plhs[9] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
   arr(double) temp_f = getArrayInput<double>(plhs[9]);
   temp_f[0] = U_calc_energy();
   //mexPrintf("%f\n", temp_f[0]);
}


//////////////////////////////////////////////////////////////////////
// 4. Utility Funcitons
//////////////////////////////////////////////////////////////////////
void IMG::U_initialize()
{
	// populate the linked lists and pointers
	for (int i=0; i<N; i++)
   {
      if (label[i]>=0)
      {
         SP_arr[label[i]]->add_pixel_init(data, i, U_check_border_pix(i), pixel_ptr[i], border_ptr[i], boundary_mask[i]);
         SP_arr[label[i]]->update_neighbors_add_self(label, i, xdim, ydim);
      }
   }

   for (int k=0; k<K; k++)
      SP_arr[k]->calculate_log_probs();
}


double IMG::U_calc_energy()
{
   double E = 0;
   int numEmpty = 0;
   for (int k=0; k<N; k++)
   {
      if (SP_arr[k]!=NULL)
      {
         SP_arr[k]->calculate_log_probs();
         E += SP_arr[k]->get_log_likelihood() + U_calc_model_order(SP_arr[k]->get_N(), SP_old[k]);
      }
      else
         numEmpty++;
   }
   E += numEmpty*(SP_new->get_log_likelihood() + U_calc_model_order(0, false));
   return E;
}

double IMG::U_calc_model_order(int size, bool is_old)
{
   if (is_old)
   {
      if (size==0)
         return 0;
      else
         return log_beta - 0.5*1.837877066409 - 0.5*log_area_var - 0.5*pow(area-size,2)/area_var;
   }
   else
   {
      if (size==0)
         return 0;
      else
         return log_alpha - 0.5*1.837877066409 - 0.5*log_area_var - 0.5*pow(area-size,2)/area_var;
   }
}


bool IMG::U_check_border_pix(int index)
{
	int x = index%xdim;
	int y = index/xdim;
	int cur_label = label[index];
	bool border = false;
	border = border || (x>0 && label[index-1]!=cur_label);
	border = border || (y>0 && label[index-xdim]!=cur_label);
	border = border || (x<xdim-1 && label[index+1]!=cur_label);
	border = border || (y<ydim-1 && label[index+xdim]!=cur_label);
	return border;
}
bool IMG::U_check_border_pix(int index, int cur_label)
{
	int x = index%xdim;
	int y = index/xdim;
	bool border = false;
	border = border || (x>0 && label[index-1]!=cur_label);
	border = border || (y>0 && label[index-xdim]!=cur_label);
	border = border || (x<xdim-1 && label[index+1]!=cur_label);
	border = border || (y<ydim-1 && label[index+xdim]!=cur_label);
	return border;
}
// --------------------------------------------------------------------------
// -- U_update_border_changed
// --   Updates all border linked lists and pointers for the neighbors of a
// -- changed pixel at index.
// --
// --   parameters:
// --     - index : the index of the recently changed pixel
// --------------------------------------------------------------------------
void IMG::U_update_border_changed(int index)
{
	int x = index%xdim;
	int y = index/xdim;
	if (x>0) U_update_border_changed_pixel(index-1);
	if (y>0) U_update_border_changed_pixel(index-xdim);
	if (x<xdim-1) U_update_border_changed_pixel(index+1);
	if (y<ydim-1) U_update_border_changed_pixel(index+xdim);
}
// --------------------------------------------------------------------------
// -- U_update_border_changed_pixel
// --   Updates all border linked lists and pointers at index
// --
// --   parameters:
// --     - index : the index of the recently changed pixel
// --------------------------------------------------------------------------
void IMG::U_update_border_changed_pixel(int index)
{
   if (label[index]>=0)
   {
      if (border_ptr[index]!=NULL && !U_check_border_pix(index))
      {
         SP_arr[label[index]]->borders.deleteNode(border_ptr[index]);
         border_ptr[index] = NULL;
      }
      else if (border_ptr[index]==NULL && U_check_border_pix(index))
         border_ptr[index] = SP_arr[label[index]]->borders.addNodeEnd(index);
   }
}

// --------------------------------------------------------------------------
// -- U_update_neighbor_list
// --   Updates super pixel neighbors used in merge
// --
// --   parameters:
// --     - neighbors : a boolean array indicating which labels are added
// --     - neighborsLL : a linked list of the neighbor indices
// --     - index : the index of the point to consider adding
// --------------------------------------------------------------------------
void IMG::U_update_neighbor_list(arr(bool) neighbors, linkedList<int> &neighborsLL, int index)
{
   int k = label[index];
	if (k>=0 && !neighbors[k])
	{
		neighborsLL.addNodeEnd(k);
		neighbors[k] = true;
	}
}

void IMG::U_relabel_SP(bool final)
{
   int last_k = 0;
   for (int k=0; k<K; k++)
   {
      if (SP_arr[k]!=NULL && (!final || SP_arr[k]->get_N()>0))
      {
         if (k!=last_k)
         {
            // move the super pixel to the empty one
            SP_arr[last_k] = SP_arr[k];
            SP_arr[k] = NULL;

            // relabel the pixels
            linkedListNode<int>* node = SP_arr[last_k]->pixels.getFirst();
            while (node != NULL)
            {
               label[node->getData()] = last_k;
               node = node->getNext();
            }
         }
         last_k++;
      }
   }
   K = last_k;
}


void IMG::U_fix_neighbors_self(int k)
{
   SP_arr[k]->neighbors.clear();
   linkedListNode<int>* node = SP_arr[k]->borders.getFirst();
   while (node!=NULL)
   {
      SP_arr[k]->update_neighbors_add_self(label, node->getData(), xdim, ydim);
      node = node->getNext();
   }
}
void IMG::U_print_neighbors(int k)
{
   linkedListNode< std::pair<int, int> >* neighborNode;
   neighborNode = SP_arr[k]->neighbors.getFirst();
   mexPrintf("%d 's NB (%d): ",k,SP_arr[k]->neighbors.getLength());
   while (neighborNode!=NULL)
   {
      mexPrintf("%d   ",neighborNode->getData().first);
      neighborNode = neighborNode->getNext();
   }
   mexPrintf("\n");
}
void IMG::U_fix_neighbors_neighbors(int k)
{
   linkedListNode< std::pair<int, int> >* neighborNode;
   neighborNode = SP_arr[k]->neighbors.getFirst();
   while (neighborNode!=NULL)
   {
      int k2 = neighborNode->getData().first;
      U_fix_neighbors_self(k2);
      neighborNode = neighborNode->getNext();
   }
}
void IMG::U_fix_neighbors_neighbors(int k, int kignore)
{
   linkedListNode< std::pair<int, int> >* neighborNode;
   neighborNode = SP_arr[k]->neighbors.getFirst();
   while (neighborNode!=NULL)
   {
      int k2 = neighborNode->getData().first;
      if (k2!=kignore)
         U_fix_neighbors_self(k2);
      neighborNode = neighborNode->getNext();
   }
}

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
void IMG::U_update_neighbors_rem(int old_label, int index)
{
   if (old_label>=0)
   {
      int x = index%xdim;
      int y = index/xdim;

      // don't do all neighbors, only the unique ones
      int llabel = (x>0 ? label[index-1] : -1);
      int ulabel = (y>0 ? label[index-xdim] : -1);
      int rlabel = (x<xdim-1 ? label[index+1] : -1);
      int dlabel = (y<ydim-1 ? label[index+xdim] : -1);

      if (llabel>=0 && llabel!=old_label)
         SP_arr[old_label]->update_neighbors_label_rem(llabel);
      if (ulabel>=0 && ulabel!=old_label && ulabel!=llabel)
         SP_arr[old_label]->update_neighbors_label_rem(ulabel);
      if (rlabel>=0 && rlabel!=old_label && rlabel!=llabel && rlabel!=ulabel)
         SP_arr[old_label]->update_neighbors_label_rem(rlabel);
      if (dlabel>=0 && dlabel!=old_label && dlabel!=llabel && dlabel!=ulabel && dlabel!=rlabel)
         SP_arr[old_label]->update_neighbors_label_rem(dlabel);

      // update the neighbors' neighbors list. (not a typo!)
      if (llabel>=0 && llabel!=old_label)
         SP_arr[llabel]->update_neighbors_label_rem_check(label, index-1, xdim, ydim, old_label);
      if (ulabel>=0 && ulabel!=old_label)
         SP_arr[ulabel]->update_neighbors_label_rem_check(label, index-xdim, xdim, ydim, old_label);
      if (rlabel>=0 && rlabel!=old_label)
         SP_arr[rlabel]->update_neighbors_label_rem_check(label, index+1, xdim, ydim, old_label);
      if (dlabel>=0 && dlabel!=old_label)
         SP_arr[dlabel]->update_neighbors_label_rem_check(label, index+xdim, xdim, ydim, old_label);
   }
}
// --------------------------------------------------------------------------
// -- U_update_neighbors_add
// --   Updates all neighbor lists when a pixel is "added". A previous
// -- call to update_neighbors_rem should be completed right before this one.
// -- The neighboring label should be changed before calling this function.
// --
// --   parameters:
// --     - index : the index bordering the removed pixel
// --------------------------------------------------------------------------
void IMG::U_update_neighbors_add(int index)
{
   int cur_label = label[index];
   if (cur_label>=0)
   {
      int x = index%xdim;
      int y = index/xdim;

      // don't do all neighbors, only the unique ones
      int llabel = (x>0 ? label[index-1] : -1);
      int ulabel = (y>0 ? label[index-xdim] : -1);
      int rlabel = (x<xdim-1 ? label[index+1] : -1);
      int dlabel = (y<ydim-1 ? label[index+xdim] : -1);
      if (llabel>=0 && llabel!=cur_label)
         SP_arr[cur_label]->update_neighbors_label_add(llabel);
      if (ulabel>=0 && ulabel!=cur_label && ulabel!=llabel)
         SP_arr[cur_label]->update_neighbors_label_add(ulabel);
      if (rlabel>=0 && rlabel!=cur_label && rlabel!=llabel && rlabel!=ulabel)
         SP_arr[cur_label]->update_neighbors_label_add(rlabel);
      if (dlabel>=0 && dlabel!=cur_label && dlabel!=llabel && dlabel!=ulabel && dlabel!=rlabel)
         SP_arr[cur_label]->update_neighbors_label_add(dlabel);

      // update the neighbors' neighbors list. (not a typo!)
      if (llabel>=0 && llabel!=cur_label)
         SP_arr[llabel]->update_neighbors_label_add_check(label, index-1, xdim, ydim, cur_label);
      if (ulabel>=0 && ulabel!=cur_label)
         SP_arr[ulabel]->update_neighbors_label_add_check(label, index-xdim, xdim, ydim, cur_label);
      if (rlabel>=0 && rlabel!=cur_label)
         SP_arr[rlabel]->update_neighbors_label_add_check(label, index+1, xdim, ydim, cur_label);
      if (dlabel>=0 && dlabel!=cur_label)
         SP_arr[dlabel]->update_neighbors_label_add_check(label, index+xdim, xdim, ydim, cur_label);
   }
}

// --------------------------------------------------------------------------
// -- U_update_neighbors_merge
// --   Updates the neighbor lists for merging two super pixels.
// -- All labels and SP_arr things should be updated *before* calling this.
// --
// --   parameters:
// --     - index : the index bordering the removed pixel
// --------------------------------------------------------------------------
void IMG::U_update_neighbors_merge(int new_label, int old_label)
{
   linkedListNode< std::pair<int, int> > *node, *node2;

   // first update the merged neighbor list
   node = SP_arr[old_label]->neighbors.getFirst();
   while (node!=NULL)
   {
      std::pair<int, int> old = node->getData();
      int neighbor_k = old.first;

      if (neighbor_k != new_label)
      {
         node2 = SP_arr[new_label]->neighbors.getFirst();
         bool found = false;
         while (node2!=NULL)
         {
            if (node2->getData().first == neighbor_k)
            {
               found = true;
               break;
            }
            node2 = node2->getNext();
         }
         if (found) // the neighbor already exists, just aggregate counts
            node2->applyFunction( &increment_neighbor_count, old );
         else // the neighbor doesn't exist.  create a new one
            SP_arr[new_label]->neighbors.addNodeEnd(old);
         node = node->getNext();
      }
   }

   // now update all neighboring neighbor lists
   node = SP_arr[new_label]->neighbors.getFirst();
   while (node!=NULL)
   {
      int neighbor_k = node->getData().first;
      SP_arr[neighbor_k]->update_neighbors_self(label, xdim, ydim);
      node = node->getNext();
   }
}

// --------------------------------------------------------------------------
// -- U_update_neighbors_split
// --   Updates the neighbor lists for merging two super pixels.
// -- All labels and SP_arr things should be updated *before* calling this.
// --
// --   parameters:
// --     - index : the index bordering the removed pixel
// --------------------------------------------------------------------------
void IMG::U_update_neighbors_split(int label1, int label2)
{
   SP_arr[label1]->update_neighbors_self(label, xdim, ydim);
   SP_arr[label2]->update_neighbors_self(label, xdim, ydim);

   linkedList<int> neighbors_done;

   linkedListNode< std::pair<int, int> > *node;
   node = SP_arr[label1]->neighbors.getFirst();
   while (node!=NULL)
   {
      int neighbor_k = node->getData().first;
      if (neighbor_k != label2)
      {
         SP_arr[neighbor_k]->update_neighbors_self(label,xdim,ydim);
         neighbors_done.addNodeEnd(neighbor_k);
      }
      node = node->getNext();
   }

   node = SP_arr[label2]->neighbors.getFirst();
   while (node!=NULL)
   {
      int neighbor_k = node->getData().first;
      if (neighbor_k != label1)
      {
         linkedListNode<int>* node2 = neighbors_done.getFirst();
         bool done = false;
         while (node2!=NULL)
         {
            if (node2->getData() == neighbor_k)
            {
               done = true;
               break;
            }
            node2 = node2->getNext();
         }
         if (!done)
            SP_arr[neighbor_k]->update_neighbors_self(label,xdim,ydim);
      }
   }
}




//////////////////////////////////////////////////////////////////////
// 6. Moving Functions
//////////////////////////////////////////////////////////////////////

void IMG::Flow_QP2()
{
   if (prev_label!=NULL)
   {
      //U_find_prev_border();
      //given label z, globally update flow variable of SP using GP approximation of Optical Flow like L2 penalty
      //global variable
      int count, tmp_ind,shift=prev_K+1;
      pair<int, int> tmp_pair;
      double lambda = 1E-3;
      arr(double) tmp_mean;
      linkedListNode< pair<int, int> >* node_p;
      linkedListNode<int>* node;
      linkedList<int> nb;

      // get the alive stuff
      num_alive = 0;
      num_dead = 0;
      for (int k=0; k<prev_K; k++)
      {
         if (!SP_arr[k]->isempty())
         {
            all2alive[k] = num_alive;
            alive2all[num_alive] = k;
            num_alive++;
         }
         else
         {
            all2dead[k] = num_dead;
            dead2all[num_dead] = k;
            num_dead++;
         }
      }

      // for computational efficiency
      // 1. Covariance matrix
      if (alive_dead_changed==true || SxySyy==NULL)
      {
         if (SxySyy!=NULL)
         {
            deallocate_memory(Sxy);
            deallocate_memory(Syy);
            deallocate_memory(SxySyy);
            deallocate_memory(obs_u);
            deallocate_memory(obs_v);
         }
         Syy = allocate_memory<double>(num_alive*num_alive);
         Sxy = allocate_memory<double>(num_alive*prev_K);
         SxySyy = allocate_memory<double>(num_alive*prev_K);
         obs_u = allocate_memory<double>(num_alive);
         obs_v = allocate_memory<double>(num_alive);

         find_SxySyy(all2alive, alive2all, all2dead, dead2all);
         alive_dead_changed = false;
      }
      /*arr(double) temp_Delta = SP_arr[0]->get_Delta_pos();
      double noise_var = temp_Delta[0];
      if (alive_dead_changed==true || SxySyy==NULL)
      {
         if (SxySyy!=NULL)
         {
            deallocate_memory(SxySyy);
            deallocate_memory(obs_u);
            deallocate_memory(obs_v);
         }
         arr(double) Syy = allocate_memory<double>(num_alive*num_alive);
         arr(double) Sxy = allocate_memory<double>(num_alive*prev_K);
         SxySyy = allocate_memory<double>(num_alive*prev_K);
         obs_u = allocate_memory<double>(num_alive);
         obs_v = allocate_memory<double>(num_alive);

         for (int k=0; k<prev_K; k++) for (int alive_k=0; alive_k<num_alive; alive_k++)
            Sxy[k*num_alive + alive_k] = prev_covariance[k*prev_K + alive2all[alive_k]];
         for (int k1=0; k1<num_alive; k1++)
         {
            for (int k2=0; k2<num_alive; k2++)
               Syy[k1*num_alive + k2] = prev_covariance[alive2all[k1]*prev_K + alive2all[k2]];
            Syy[k1*num_alive + k1] += noise_var;
         }

         // cholesky stored in Syy
         gsl_matrix_view Syy_view = gsl_matrix_view_array(Syy, num_alive, num_alive);
         gsl_linalg_cholesky_decomp(&Syy_view.matrix);
         gsl_linalg_cholesky_invert(&Syy_view.matrix);

         gsl_matrix_view Sxy_view = gsl_matrix_view_array(Sxy, prev_K, num_alive);
         gsl_matrix_view SxySyy_view = gsl_matrix_view_array(SxySyy, prev_K, num_alive);

         // Sxy*Syy_chol
         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Sxy_view.matrix, &Syy_view.matrix, 0.0, &SxySyy_view.matrix);

         deallocate_memory(Syy);
         deallocate_memory(Sxy);
         alive_dead_changed = false;
      }*/

      // populate the observation
      int all_k1;
      for (int k1=0; k1<num_alive; k1++)
      {
         all_k1 = alive2all[k1];

         arr(double) new_pos_mean = SP_arr[all_k1]->get_mean_pos();
         arr(double) prev_velocity = SP_arr[all_k1]->get_prev_v();
         obs_u[k1] = new_pos_mean[0] - prev_pos_mean[all_k1*2] - prev_velocity[0];
         obs_v[k1] = new_pos_mean[1] - prev_pos_mean[all_k1*2+1] - prev_velocity[1];
         //mexPrintf("Diff: %d, %f , %f\n", k1, obs_u[k1], obs_v[k1]);
      }

      // now use GSL to compute Sxy*(Syy)^-1 * obs_X
      // (SxySyy) * obs_X
      /*gsl_matrix_view SxySyy_view = gsl_matrix_view_array(SxySyy, prev_K, num_alive);
      gsl_blas_dgemv (CblasNoTrans, 1.0, &SxySyy_view.matrix, &obs_u_view.vector, 0.0, ugsl);
      gsl_blas_dgemv (CblasNoTrans, 1.0, &SxySyy_view.matrix, &obs_v_view.vector, 0.0, vgsl);*/

      arr(double) temp_vec = allocate_memory<double>(num_alive);
      gsl_matrix_view Sxy_view = gsl_matrix_view_array(Sxy, prev_K, num_alive);
      gsl_matrix_view Syy_view = gsl_matrix_view_array(Syy, num_alive, num_alive);
      gsl_vector_view obs_u_view = gsl_vector_view_array(obs_u, num_alive);
      gsl_vector_view obs_v_view = gsl_vector_view_array(obs_v, num_alive);
      gsl_vector_view temp_vec_view = gsl_vector_view_array(temp_vec, num_alive);

      gsl_blas_dgemv (CblasNoTrans, 1.0, &Syy_view.matrix, &obs_u_view.vector, 0.0, &temp_vec_view.vector);
      gsl_blas_dgemv (CblasNoTrans, 1.0, &Sxy_view.matrix, &temp_vec_view.vector, 0.0, ugsl);
      gsl_blas_dgemv (CblasNoTrans, 1.0, &Syy_view.matrix, &obs_v_view.vector, 0.0, &temp_vec_view.vector);
      gsl_blas_dgemv (CblasNoTrans, 1.0, &Sxy_view.matrix, &temp_vec_view.vector, 0.0, vgsl);

      deallocate_memory(temp_vec);

      // copy over the new parameters
      arr(double) temp_f = allocate_memory<double>(2);
      for (int i=0;i<prev_K;i++)
      {
         if (SP_arr[i]!=NULL)
         {
            temp_f[0] = gsl_vector_get(ugsl,i);
            temp_f[1] = gsl_vector_get(vgsl,i);
            arr(double) prev_velocity = SP_arr[i]->get_prev_v();
            temp_f[0] += prev_velocity[0];
            temp_f[1] += prev_velocity[1];
            SP_arr[i]->set_flow(temp_f);
            //SP_arr[i]->set_meansum_pos(temp_f, prev_pos_mean+i*2);

            //mexPrintf("flow: %d, %f , %f \n",i, gsl_vector_get(ugsl,i),gsl_vector_get(vgsl,i));
            //SP_arr[i]->set_flow(gsl_vector_get(u,i), gsl_vector_get(v,i));
         }
         SP_changed[i] = true;
      }

      // free all the memory
      deallocate_memory(temp_f);
   }
}

void IMG::find_SxySyy(arr(int) all2alive, arr(int) alive2all, arr(int) all2dead, arr(int) dead2all)
{
   //temp_Syy is prev_K x prev_K
   // dead2all
   // temp_d is prev_K x 1
   // still_alive is prev_K x 1
   // assume memory allocated correctly for SxySyy
   memset(temp_still_alive, true, sizeof(bool)*prev_K);
   memcpy(temp_Syy, prev_precision, sizeof(double)*prev_K*prev_K);
   for (int k=0; k<num_dead; k++)
   {
      int dead_k = dead2all[k];
      double c = temp_Syy[dead_k + dead_k*prev_K];
      temp_still_alive[dead_k] = false;
      /*for (int i=0; i<prev_K; i++)
         if (temp_still_alive[i])
            temp_d[i] = temp_Syy[dead_k + i*prev_K];*/
      //only populate half of the matrix
      for (int i=0; i<dead_k; i++)
         if (temp_still_alive[i])
            temp_d[i] = temp_Syy[dead_k + i*prev_K];
      int offset = dead_k*prev_K;
      for (int i=dead_k+1; i<prev_K; i++)
         if (temp_still_alive[i])
            temp_d[i] = temp_Syy[i + offset];


      for (int i=0; i<prev_K; i++) if (temp_still_alive[i])
      {
         int offset = i*prev_K;
         double temp_di = temp_d[i];
         temp_Syy[i + offset] -= (pow(temp_di,2)/c);
         temp_di /= c;
         for (int j=i+1; j<prev_K; j++) if (temp_still_alive[j])
            temp_Syy[j + offset] -= (temp_di*temp_d[j]);
      }
   }


   for (int k=0; k<prev_K; k++) for (int alive_k=0; alive_k<num_alive; alive_k++)
      Sxy[k*num_alive + alive_k] = prev_covariance[k*prev_K + alive2all[alive_k]];
   for (int k1=0; k1<num_alive; k1++)
   {
      int alive_k1 = alive2all[k1];
      int offset = alive_k1*prev_K;
      int offset_alive = k1*num_alive;
      Syy[k1 + offset_alive] = temp_Syy[alive_k1 + offset];
      for (int k2=k1+1; k2<num_alive; k2++)
      {
         int alive_k2 = alive2all[k2];
         Syy[k1 + k2*num_alive] = temp_Syy[alive_k2 + offset];
         Syy[k2 + offset_alive] = temp_Syy[alive_k2 + offset];
      }
   }

   gsl_matrix_view Syy_view = gsl_matrix_view_array(Syy, num_alive, num_alive);
   gsl_matrix_view Sxy_view = gsl_matrix_view_array(Sxy, prev_K, num_alive);
   gsl_matrix_view SxySyy_view = gsl_matrix_view_array(SxySyy, prev_K, num_alive);

   // Sxy*Syy^-1
   //gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Sxy_view.matrix, &Syy_view.matrix, 0.0, &SxySyy_view.matrix);
   //gsl_blas_dsymm(CblasRight, CblasUpper, 1.0, &Syy_view.matrix, &Sxy_view.matrix, 0.0, &SxySyy_view.matrix);
}


// --------------------------------------------------------------------------
// -- link_cost
// --   calculates the energy of linking SP_arr[ko] (old) to SP_arr[kn] (new)
// --
// --   parameters:
// --     - ko : the old k
// --     - kn : the new k
// --------------------------------------------------------------------------
double IMG::link_cost(int ko, int kn)
{
   double logprob;
   if (ko<prev_K && kn<K) // matched SPs
   {
      logprob = SP_arr[kn]->log_likelihood_switch_app_prior(SP_arr[ko]);
      logprob += U_calc_model_order(SP_arr[kn]->get_N(), true);
   }
   else if (ko>=prev_K) // new SP
   {
      logprob = SP_arr[kn]->log_likelihood_switch_app_prior(SP_new);
      logprob += U_calc_model_order(SP_arr[kn]->get_N(), false);
   }
   else if (kn>=K) // dead SP
   {
      logprob = SP_new->log_likelihood_switch_app_prior(SP_arr[ko]);
   }
   else
   {
      mexErrMsgTxt("Can't link a dead SP to a new SP!");
   }
   return logprob;
}

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
double IMG::move_local_calc_delta(int index, int new_k, bool add, double& max_prob, int& max_k)
{
   double prob = 0;
   if (new_k<0)
   {
      prob = boundary_mask[index] ? -mxGetInf() : dummy_log_prob;
   }
   else
   {
      SP* temp;
      bool is_old;
      if (SP_arr[new_k]==NULL) // new super pixel
      {
         temp = SP_new;
         is_old = false;
      }
      else
      {
         temp = SP_arr[new_k];
         is_old = SP_old[new_k];
      }

      int N = temp->get_N();
      if (add)
      {
         prob += temp->log_likelihood_test_point(data + index*5, boundary_mask[index]);
         prob -= temp->get_log_likelihood();

         prob += U_calc_model_order(N+1, is_old);
         prob -= U_calc_model_order(N, is_old);
      }
      else
      {
         prob += temp->get_log_likelihood();
         prob -= temp->log_likelihood_test_point_rem(data + index*5, boundary_mask[index]);

         prob += U_calc_model_order(N, is_old);
         prob -= U_calc_model_order(N-1, is_old);
      }
   }

	if (prob>max_prob+1e-10 || max_k==-2)
	{
		max_prob = prob;
		max_k = new_k;
	}
}



// --------------------------------------------------------------------------
// -- move_local_calc_neigbor
// --   calculates the probability of assigning the pixel at index to the
// -- cluster of nindex (neighbor index). Updates the max_prob and max_k if
// -- the new_app value is greater than the old one
// --
// --   parameters:
// --     - index : the new point to add
// --     - nindex : the neighbor to add this point to
// --     - max_prob : the maximum probability of all neighbors
// --     - max_k : the super pixel index of the maximum probability
// --------------------------------------------------------------------------
double IMG::move_local_calc_delta_MM(int index, int new_k, double& max_prob, int& max_k)
{
   double prob;
   prob = SP_arr[new_k]->log_likelihood_test_point_MM(data + index*5);

   if (prob>max_prob+1e-10 || max_k<0)
   {
      max_prob = prob;
      max_k = new_k;
   }
}




double IMG::move_switch_calc_delta(SP* &oldSP, SP* &newSP)
{
   double logprob = oldSP->log_likelihood_switch_prior(newSP);
   logprob -= oldSP->get_log_likelihood();

   bool old_isold = oldSP->is_old;
   bool new_isold = newSP->is_old;
   double Ni = oldSP->get_N();

   logprob += U_calc_model_order(Ni, new_isold);
   logprob -= U_calc_model_order(Ni, old_isold);

   return logprob;
}



bool IMG::move_switch_IMG()
{
   linkedList<int> empty_SPs;
   for (int k=0; k<K; k++)
      if (SP_arr[k]->isempty() && SP_old[k])
         empty_SPs.addNodeEnd(k);

   bool changed = false;
   for (int k=0; k<K; k++)
   {
      if (!(SP_arr[k]->isempty()))
      {
         // if old, check to see a new or unused old
         // if new, check to see if any unused old
         int best_k = -1;
         double best_energy = 0;
         double delta;
         if (SP_old[k])
         {
            delta = move_switch_calc_delta(SP_arr[k], SP_new);
            if (delta > best_energy)
            {
               best_k = K;
               best_energy = delta;
            }
         }

         linkedListNode<int>* node = empty_SPs.getFirst();
         linkedListNode<int>* best_node = NULL;
         while (node!=NULL)
         {
            int test_k = node->getData();
            delta = move_switch_calc_delta(SP_arr[k], SP_arr[test_k]);
            if (delta > best_energy)
            {
               best_k = test_k;
               best_energy = delta;
               best_node = node;
            }
            node = node->getNext();
         }

         // switch with best label
         if (best_k>=0 && best_energy>0)
         {
            // change the labels
            node = SP_arr[k]->pixels.getFirst();
            while (node != NULL)
            {
               int index = node->getData();
               label[index] = best_k;
               node = node->getNext();
            }
            SP_changed[k] = true;
            SP_changed[best_k] = true;

            // make room for the new one
            if (best_k==K)
            {
               if (K>=N)
                  mexErrMsgTxt("Ran out of space!");
               if (SP_arr[K]==NULL)
                  SP_arr[K] = new SP(new_pos, new_app, max_UID++);
               K++;
            }
            else
               alive_dead_changed = true;
            SP_arr[best_k]->merge_with(SP_arr[k], label, border_ptr, xdim, ydim);

            // delete it if it was a new one
            if (!SP_old[k])
            {
               delete SP_arr[k];
               SP_arr[k] = NULL;
            }
            else
            {
               // add it to the search list
               empty_SPs.addNodeEnd(k);
               alive_dead_changed = true;
            }

            // remove it from the search list
            if (SP_old[best_k])
               empty_SPs.deleteNode(best_node);


            // now update the neighbors list
            U_fix_neighbors_self(best_k);
            // update the neighbors' neighbors
            U_fix_neighbors_neighbors(best_k, N+1);
         }
      }
   }
}




// --------------------------------------------------------------------------
// -- move_local
// --   finds the optimal local joint move of labels and parameters. Chooses
// -- a super pixel at random and loops through and updates its borders.
// --------------------------------------------------------------------------
bool IMG::move_local_IMG()
{
   bool changed = false;

   // temporary neighborhood
   bool neighborhood[9];

   // choose a random order of super pixels
   const size_t Nsp = K;
   const gsl_rng_type *T;
   gsl_rng *r;
   gsl_rng_env_setup();
   gsl_permutation *perm = gsl_permutation_alloc(Nsp);
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   gsl_rng_set(r, (unsigned long)(rand()));
   gsl_permutation_init(perm);
   gsl_ran_shuffle(r, perm->data, Nsp, sizeof(size_t));

   if (K>N)
      mexErrMsgTxt("Ran out of space!\n");

   // MM STEP
   /*for (int k=0; k<K; k++)
   {
      SP_arr[k]->get_mean_pos(true);
      SP_arr[k]->get_mean_app();
   }*/

   for (int ki=0; ki<Nsp; ki++)
   {
      int k = perm->data[ki];

      // find a nonempty super pixel
      if (SP_arr[k]!=NULL && !(SP_arr[k]->isempty()) && SP_changed[k])
      {
         SP_changed[k] = false;
         // loop through borders
         int length = SP_arr[k]->borders.getLength();
         int i = 0;

         linkedList<int> dummy_indices;
         while (i<length*5)
         {
            if (k>=0 && (SP_arr[k]==NULL || SP_arr[k]->isempty()))
               break;

            int index;
            if (i<length)
               index = SP_arr[k]->borders.popFirst();
            else
            {
               if (dummy_indices.isempty())
                  break;
               else
                  index = dummy_indices.popFirst();
               k = label[index];
               if (k>=0)
               {
                  i++;
                  continue;
               }
            }
            i++;

            // check the topology for this pixel
            if (!check_topology(index, label, neighborhood, xdim, ydim, T4Table))
            {
               if (k>=0)
                  border_ptr[index] = SP_arr[k]->borders.addNodeEnd(index);
            }
            else
            {
               // temporarily remove this data point from the SP
               double max_prob = -mxGetInf();
               int max_k = -2;
               int x = index%xdim;
               int y = index/xdim;

               // the current k has to be a possible choice
               move_local_calc_delta(index, label[index], false, max_prob, max_k);
               // a new k is also a possible choice
               //move_local_calc_delta(index, K, true, max_prob, max_k);

               if (!boundary_mask[index])
                  move_local_calc_delta(index, -1, true, max_prob, max_k);

               // find which k's we can move to
               if (x>0 && label[index-1]!=k) move_local_calc_delta(index, label[index-1], true, max_prob, max_k);
               if (y>0 && label[index-xdim]!=k) move_local_calc_delta(index, label[index-xdim], true, max_prob, max_k);
               if (x<xdim-1 && label[index+1]!=k) move_local_calc_delta(index, label[index+1], true, max_prob, max_k);
               if (y<ydim-1 && label[index+xdim]!=k) move_local_calc_delta(index, label[index+xdim], true, max_prob, max_k);

               /*move_local_calc_delta_MM(index, label[index], max_prob, max_k);
               if (x>0 && label[index-1]!=k) move_local_calc_delta_MM(index, label[index-1], max_prob, max_k);
               if (y>0 && label[index-xdim]!=k) move_local_calc_delta_MM(index, label[index-xdim], max_prob, max_k);
               if (x<xdim-1 && label[index+1]!=k) move_local_calc_delta_MM(index, label[index+1], max_prob, max_k);
               if (y<ydim-1 && label[index+xdim]!=k) move_local_calc_delta_MM(index, label[index+xdim], max_prob, max_k);*/

               if (max_k!=k)
               {
                  changed = true;
                  // update the labels... it moves from k->max_k
                  label[index] = max_k;
                  if (max_k>=K)
                  {
                     // creating a new one
                     if (K>=N-1)
                        mexErrMsgTxt("Ran out of space!\n");
                     SP_arr[K] = new SP(new_pos, new_app, max_UID++);
                     K++;
                  }

                  // update the neighbors lists
                  U_update_neighbors_rem(k, index);
                  U_update_neighbors_add(index);

                  // set the correct SP_changed variables of all neighbors
                  if (k<0)
                  {
                     SP_changed[max_k] = true;
                  }
                  else
                  {
                     linkedListNode< pair<int, int> >* node = SP_arr[k]->neighbors.getFirst();
                     while (node != NULL)
                     {
                        int neighbor_k = node->getData().first;
                        SP_changed[neighbor_k] = true;
                        node = node->getNext();
                     }
                     if (max_k>=0)
                     {
                        node = SP_arr[max_k]->neighbors.getFirst();
                        while (node != NULL)
                        {
                           int neighbor_k = node->getData().first;
                           SP_changed[neighbor_k] = true;
                           node = node->getNext();
                        }
                     }
                  }


                  // update all border lists for neighbors
                  U_update_border_changed(index);
                  if (k>=0)
                     SP_arr[k]->rem_pixel(data, index, pixel_ptr[index], boundary_mask[index]);

                  if (k>=0 && SP_arr[k]->isempty())
                  {
                     if (SP_old[k])
                        alive_dead_changed = true;
                     else
                     {
                        delete SP_arr[k];
                        SP_arr[k] = NULL;
                        SP_changed[k] = false;
                     }
                  }

                  // add this point to the maximum SP
                  if (max_k>=0)
                     SP_arr[max_k]->add_pixel(data, index, U_check_border_pix(index), pixel_ptr[index], border_ptr[index], boundary_mask[index]);
                  else
                  {
                     pixel_ptr[index] = NULL;
                     border_ptr[index] = NULL;
                  }
               }
               else
               {
                  if (k>=0)
                     border_ptr[index] = SP_arr[k]->borders.addNodeEnd(index);
               }
               if (x>0 && label[index-1]<0)
                  dummy_indices.addNodeEnd(index-1);
               if (y>0 && label[index-xdim]<0)
                  dummy_indices.addNodeEnd(index-xdim);
               if (x<xdim-1 && label[index+1]<0)
                  dummy_indices.addNodeEnd(index+1);
               if (y<ydim-1 && label[index+xdim]<0)
                  dummy_indices.addNodeEnd(index+xdim);
            }
         }
      }
   }

   if (!changed)
   {
      for (int k=0; k<K; k++)
         SP_changed[k] = false;
   }

   gsl_permutation_free(perm);
   gsl_rng_free(r);
   return changed;
}


// --------------------------------------------------------------------------
// -- move_merge_calc_delta
// --   calculates the probability of assigning the pixel at index to the
// -- cluster of nindex (neighbor index).
// --
// --   parameters:
// --     - index : the new point to add
// --     - nindex : the neighbor to add this point to
// --------------------------------------------------------------------------
double IMG::move_merge_calc_delta(int k, int merge_k)
{
	double prob = SP_arr[k]->log_likelihood_test_merge(SP_arr[merge_k]) + SP_arr[merge_k]->get_log_likelihood_empty();
	prob -= (SP_arr[k]->get_log_likelihood() + SP_arr[merge_k]->get_log_likelihood());

   prob += U_calc_model_order(SP_arr[k]->get_N() + SP_arr[merge_k]->get_N(), SP_old[k]);
   prob += U_calc_model_order(0, SP_old[merge_k]);
   prob -= U_calc_model_order(SP_arr[k]->get_N(), SP_old[k]);
   prob -= U_calc_model_order(SP_arr[merge_k]->get_N(), SP_old[merge_k]);

	//prob += calc_log_label(SP_arr[k]->get_N() + SP_arr[merge_k]->get_N(), alpha, epsilon);
	//prob -= (SP_arr[k]->get_log_label() + SP_arr[merge_k]->get_log_label());

   return prob;
}

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
double IMG::move_split_calc_delta(SP* SP1, SP* SP2, SP* new_SP1, SP* new_SP2, bool SP1_old, bool SP2_old)
{
   double prob = 0;
   prob += SP1->log_likelihood_test_merge(new_SP1);
   prob += SP2->log_likelihood_test_merge(new_SP2);
   prob -= SP1->log_likelihood_test_merge(new_SP1, new_SP2);
   prob -= SP2->get_log_likelihood();

   // split
   prob += U_calc_model_order(SP1->get_N() + new_SP1->get_N(), SP1_old);
   prob += U_calc_model_order(SP2->get_N() + new_SP2->get_N(), SP2_old);

   // not split
   prob -= U_calc_model_order(SP1->get_N() + new_SP1->get_N() + new_SP2->get_N(), SP1_old);
   prob -= U_calc_model_order(SP2->get_N(), SP2_old);

   return prob;
}

// --------------------------------------------------------------------------
// -- move_merge
// --   finds the optimal merge between all pairs of super pixels
// --------------------------------------------------------------------------
void IMG::move_merge_IMG()
{
   // choose a random order of super pixels
   const size_t Nsp = K;
   const gsl_rng_type *T;
   gsl_rng *r;
   gsl_rng_env_setup();
   gsl_permutation *perm = gsl_permutation_alloc(Nsp);
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   gsl_rng_set(r, (unsigned long)(rand()));
   gsl_permutation_init(perm);
   gsl_ran_shuffle(r, perm->data, Nsp, sizeof(size_t));
   

   arr(bool) neighbors = allocate_memory<bool>(K, false);
   bool merged = false;
   for (int ki=0; ki<K; ki++)
   {
      int k = perm->data[ki];

      // find a nonempty super pixel
      if (SP_arr[k]!=NULL && !(SP_arr[k]->isempty()))
      {
         linkedList<int> neighborsLL;
         // find all bordering super pixels
         linkedListNode<int>* node = SP_arr[k]->borders.getFirst();
         while (node != NULL)
         {
            int index = node->getData();
            int x = index%xdim;
            int y = index/xdim;

            if (x>0 && label[index-1]!=k) U_update_neighbor_list(neighbors, neighborsLL, index-1);
            if (y>0 && label[index-xdim]!=k) U_update_neighbor_list(neighbors, neighborsLL, index-xdim);
            if (x<xdim-1 && label[index+1]!=k) U_update_neighbor_list(neighbors, neighborsLL, index+1);
            if (y<ydim-1 && label[index+xdim]!=k) U_update_neighbor_list(neighbors, neighborsLL, index+xdim);

            node = node->getNext();
         }

         double max_E = -DBL_MAX;
         int max_k = -1;

         // loop through all neighbors
         node = neighborsLL.getFirst();
         while (node != NULL)
         {
            int merge_k = node->getData();

            // calculate the energy
            if (merge_k>=0)
            {
               double new_E = move_merge_calc_delta(k, merge_k);
               if (new_E > max_E || max_k==-1)
               {
                  max_E = new_E;
                  max_k = merge_k;
               }

               // fix neighbors for the next super pixel
               neighbors[node->getData()] = false;
            }
            node = node->getNext();
         }

         // merge if it increases energy
         if (max_E>0)
         {
            merged = true;
            // change the labels
            node = SP_arr[max_k]->pixels.getFirst();
            while (node != NULL)
            {
               int index = node->getData();
               label[index] = k;
               node = node->getNext();
            }
            SP_changed[k] = true;
            SP_changed[max_k] = true;

            SP_arr[k]->merge_with(SP_arr[max_k], label, border_ptr, xdim, ydim);
            if (!SP_old[max_k])
            {
               delete SP_arr[max_k];
               SP_arr[max_k] = NULL;
            }
            else
               alive_dead_changed = true;
         }
         /*else if (SP_arr[k]->get_N()==1)
         {
            mexPrintf("k=%d\n", k);
            mexPrintf("E=%e\n", max_E);
            mexErrMsgTxt("Not merging one pixel region");
         }*/
      }
   }

   gsl_permutation_free(perm);
   gsl_rng_free(r);
}


void IMG::move_IMG()
{
   // choose a random order of super pixels
   const size_t Nsp = K;
   const gsl_rng_type *T;
   gsl_rng *r;
   gsl_rng_env_setup();
   gsl_permutation *perm = gsl_permutation_alloc(Nsp);
   T = gsl_rng_default;
   r = gsl_rng_alloc(T);
   gsl_rng_set(r, (unsigned long)(rand()));
   gsl_permutation_init(perm);
   gsl_ran_shuffle(r, perm->data, Nsp, sizeof(size_t));
   arr(bool) neighbors = allocate_memory<bool>(K, false);

   int pre_K = K;
   int split_thres =  N/K;
   split_thres = 2;

   arr(double) energies = allocate_memory<double>(K);
   double max_energy = -mxGetInf();
   double min_energy = mxGetInf();
   double mean_area = 0;
   for (int k=0; k<pre_K; k++)
   {
      double temp = (SP_arr[k]->get_log_likelihood() + U_calc_model_order(SP_arr[k]->get_N(), SP_old[k]))/SP_arr[k]->get_N();
      energies[k] = temp;
      if (temp>max_energy) max_energy = temp;
      if (temp<min_energy) min_energy = temp;
      //mexPrintf("%d\t%e\n", k, temp);
      mean_area += SP_arr[k]->get_N();
   }
   mean_area /= pre_K;
   double threshold = min_energy + (max_energy-min_energy)*0.2;

   for (int ki=0; ki<pre_K; ki++)
   {
      int k = perm->data[ki];
      //if (false && SP_changed[k] && SP_arr[k]!=NULL && (SP_arr[k]->get_N()>split_thres))
      if (!SP_arr[k]->has_no_neighbors() && (SP_arr[k]->get_N()>mean_area || energies[k] < threshold) && SP_changed[k])
      {
         //printf("to split: %d,%d\n",ki,k);
         move_split_SP(k);
         //move_refine_SP(k);
          //if(ki>40){ break;}
      }
   }
   gsl_permutation_free(perm);
   gsl_rng_free(r);
   deallocate_memory(energies);
}



int IMG::move_local_SP(int k)
{
   // temporary neighborhood
   bool neighborhood[9];

   // choose a random order of super pixels

   // find a nonempty super pixel
   int changed = 0;
   if (SP_arr[k]!=NULL && !(SP_arr[k]->isempty()))
   {
      // loop through borders
      int length = SP_arr[k]->borders.getLength();
      for (int i=0; i<length; i++)
      {
         if (SP_arr[k]->isempty())
            break;
         int index = SP_arr[k]->borders.popFirst();

         // check the topology for this pixel
         if (false && !check_topology(index, label, neighborhood, xdim, ydim, T4Table))
            border_ptr[index] = SP_arr[k]->borders.addNodeEnd(index);
         else
         {
            // temporarily remove this data point from the SP
            // SP_arr[k]->rem_pixel(data, index, pixel_ptr[index]);

            double max_prob = -10^10;
            int max_k = -1;
            int x = index%xdim;
            int y = index/xdim;

            // the current k has to be a possible choice
            move_local_calc_delta(index, label[index], false, max_prob, max_k);
            // a new k is also a possible choice
            move_local_calc_delta(index, K, true, max_prob, max_k);

            // find which k's we can move to
            if (x>0 && label[index-1]!=k) move_local_calc_delta(index, label[index-1], true, max_prob, max_k);
            if (y>0 && label[index-xdim]!=k) move_local_calc_delta(index, label[index-xdim], true, max_prob, max_k);
            if (x<xdim-1 && label[index+1]!=k) move_local_calc_delta(index, label[index+1], true, max_prob, max_k);
            if (y<ydim-1 && label[index+xdim]!=k) move_local_calc_delta(index, label[index+xdim], true, max_prob, max_k);

            if (max_k!=k)
            {
               // update the labels
               label[index] = max_k;
               if (max_k>=K)
               {
                  SP_arr[K] = new SP(new_pos, new_app, max_UID++);
                  K++;
               }

               // update all border lists for neighbors
               U_update_border_changed(index);

               SP_arr[k]->rem_pixel(data, index, pixel_ptr[index], border_ptr[index], boundary_mask[index]);
               if (SP_arr[k]->isempty())
               {
                  if (SP_old[k])
                     alive_dead_changed = true;
                  else
                  {
                     delete SP_arr[k];
                     SP_arr[k] = NULL;
                  }
               }

               // add this point to the maximum SP
               SP_arr[max_k]->add_pixel(data, index, U_check_border_pix(index), pixel_ptr[index], border_ptr[index], boundary_mask[index]);
            }
            else
               border_ptr[index] = SP_arr[k]->borders.addNodeEnd(index);
         }
      }
   }
   return changed;
}

void IMG::U_find_border_SP(int k, arr(bool) neighbors, linkedList<int> &neighborsLL)
{
   // find the bordering pixels of one super pixel
   linkedListNode<int>* node = SP_arr[k]->borders.getFirst();
   //printf("fff:%d,%d,%d\n",k,SP_arr[k]->N,SP_arr[k]->borders.getLength());
   while (node != NULL)
   {
      int index = node->getData();
      int x = index%xdim;
      int y = index/xdim;

      if (x>0 && label[index-1]!=k) U_update_neighbor_list(neighbors, neighborsLL, index-1);
      if (y>0 && label[index-xdim]!=k) U_update_neighbor_list(neighbors, neighborsLL, index-xdim);
      if (x<xdim-1 && label[index+1]!=k) U_update_neighbor_list(neighbors, neighborsLL, index+1);
      if (y<ydim-1 && label[index+xdim]!=k) U_update_neighbor_list(neighbors, neighborsLL, index+xdim);
         // if(index+xdim<N) printf("mmii:%d,%d,%d,%d,%d\n",x,y,k,label[index+xdim],neighborsLL.getLength());
      node = node->getNext();
   }
}
void IMG::move_merge_SP_propose(int k,arr(bool) neighbors,double &max_E,int &max_k)
{
   linkedList<int> neighborsLL;
   U_find_border_SP(k,neighbors,neighborsLL);
   // loop through all neighbors
   linkedListNode<int>* node = neighborsLL.getFirst();
   double tmp_E;
   int merge_k;
   while (node != NULL)
   {
      merge_k = node->getData();
      // calculate the energy
      tmp_E = move_merge_calc_delta(k, merge_k);
         if(tmp_E>max_E){
            max_E = tmp_E;
            max_k = merge_k;
         }
      // fix neighbors for the next super pixel
      neighbors[node->getData()] = false;
      node = node->getNext();
   }
   if (max_k==-1)
      mexErrMsgTxt("find_bestnb: No neighbour found\n");
}
void IMG::move_merge_SP_propose_region(int k,arr(bool) neighbors,linkedList<int> &check_labels,double &max_E,int &max_k)
{
   linkedList<int> neighborsLL;
   U_find_border_SP(k,neighbors,neighborsLL);

   // loop through all neighbors
   linkedListNode<int>* node = neighborsLL.getFirst();
   double tmp_E = 0;
   int merge_k;
   while (node != NULL)
   {
      merge_k = node->getData();
      //printf("inside %d,%d,%d\n",merge_k,k,check_labels.getLength());
      if(U_FindArray(merge_k,check_labels)!=-1)
      {
         // calculate the energy
         tmp_E = move_merge_calc_delta(k, merge_k);
         if(tmp_E>max_E)
         {
               max_E = tmp_E;
               max_k = merge_k;
         }
      }
      // fix neighbors for the next super pixel
      neighbors[merge_k] = false;
      node = node->getNext();
   }
   if (max_k==-1)
   {
      neighborsLL.print();
      mexPrintf("k=%d\n", debug);
      mexErrMsgTxt("find_bestnb_region: No neighbour found\n");
   }
}
void IMG::move_merge_SP(int k,arr(bool) neighbors)
{
   // find a nonempty super pixel
   if (SP_arr[k]!=NULL && !(SP_arr[k]->isempty()))
   {
      double max_E = -DBL_MAX;
      int max_k = -1;
      move_merge_SP_propose( k, neighbors,max_E,max_k);
      // merge if it increases energy
      if (max_E>0)
      {
         // change the labels
         linkedListNode<int>* node = SP_arr[max_k]->pixels.getFirst();
         while (node != NULL)
         {
            int index = node->getData();
            label[index] = k;
            node = node->getNext();
         }
         SP_arr[k]->merge_with(SP_arr[max_k], label, border_ptr, xdim, ydim);
         if (SP_old[max_k])
         {
            delete SP_arr[max_k];
            SP_arr[max_k] = NULL;
         }
      }
   }
}






void IMG::move_split_SP_propose(int index, int num_SP, int option, double &max_E,int &ksplit,int* new_ks)
{
   int num_iter = 5,SP_bbox[4]={xdim,0,ydim,0};
   // 1. Kmeans++
   bool broken = U_Kmeans_plusplus(index, SP_bbox, num_SP, num_iter, option);

   if (broken)
   {
      max_E = -DBL_MAX;
      ksplit = -2;
      return;
   }

   // 2. create two new SP from the old one
   // label matrix is equivalent to SP.pixels
   // we will mess around label matrix for new proposal
   // and recover it from SP.pixels if it doesn't work out
   // merge small connected component in the label matrix
   U_connect_newSP(SP_bbox, num_SP);

   // new super pixels in K and K+1... old super pixel in index
   SP_arr[index]->empty(false);
   double E;
   if (option==-1)
   {
      // (index, new);
      E = move_split_calc_delta(SP_arr[index], SP_new, SP_arr[K], SP_arr[K+1], SP_old[index], false);
      if (E>max_E)
      {
         max_E = E;
         new_ks[0] = K;
         new_ks[1] = K+1;
         ksplit = -1;
      }
      // (new, index)
      E = move_split_calc_delta(SP_arr[index], SP_new, SP_arr[K+1], SP_arr[K], false, SP_old[index]);
      if (E>max_E)
      {
         max_E = E;
         new_ks[0] = K+1;
         new_ks[1] = K;
         ksplit = -1;
      }
      // (index, old_empty) && (old_empty, index)
      for (int ktest=0; ktest<K; ktest++) if (ktest!=index && SP_arr[ktest]!=NULL && SP_arr[ktest]->get_N()==0)
      {
         E = move_split_calc_delta(SP_arr[index], SP_arr[ktest], SP_arr[K], SP_arr[K+1], SP_old[index], SP_old[ktest]);
         if (E>max_E)
         {
            max_E = E;
            new_ks[0] = K;
            new_ks[1] = K+1;
            ksplit = ktest;
         }

         E = move_split_calc_delta(SP_arr[index], SP_arr[ktest], SP_arr[K+1], SP_arr[K], SP_old[index], SP_old[ktest]);
         if (E>max_E)
         {
            max_E = E;
            new_ks[0] = K+1;
            new_ks[1] = K;
            ksplit = ktest;
         }
      }
      /*mexPrintf("difference=%e\n", max_E + move_merge_calc_delta(K,K+1));
      max_E = -move_merge_calc_delta(K, K+1);
      new_ks[0] = K;
      new_ks[1] = K+1;
      ksplit = -1;*/

   }
   else
   {
      // for refine_move
      // split into SP_arr[index] and SP_arr[option]
      double E_1_K1;
      double E_1_K2;
      if (!SP_old[index] && !SP_old[option])
      {
         arr(double) SP1_total_pos = SP_arr[K-2]->get_total_pos();
         double N1 = SP_arr[K-2]->get_N();
         arr(double) SP2_total_pos = SP_arr[K-1]->get_total_pos();
         double N2 = SP_arr[K-1]->get_N();
         arr(double) SPK1_total_pos = SP_arr[K]->get_total_pos();
         double NK1 = SP_arr[K]->get_N();
         arr(double) SPK2_total_pos = SP_arr[K+1]->get_total_pos();
         double NK2 = SP_arr[K+1]->get_N();

         E_1_K1 = 0;
         E_1_K2 = 0;
         for (int d=0; d<2; d++)
         {
            double temp = SPK1_total_pos[d]/NK1 - SP1_total_pos[d]/N1;
            E_1_K1 += temp*temp;
            temp = SPK2_total_pos[d]/NK2 - SP2_total_pos[d]/N2;
            E_1_K1 += temp*temp;

            temp = SPK1_total_pos[d]/NK1 - SP2_total_pos[d]/N2;
            E_1_K2 += temp*temp;
            temp = SPK2_total_pos[d]/NK2 - SP1_total_pos[d]/N1;
            E_1_K2 += temp*temp;
         }

         if (E_1_K1 < E_1_K2)
         {
            //if ((index==177 && option==493) || (index==493 && option==177)){mexPrintf("1) E1K1=%e, E1K2=%e\n", E_1_K1, E_1_K2);drawnow();}
            new_ks[0] = K;
            new_ks[1] = K+1;
            max_E = move_split_calc_delta(SP_arr[index], SP_arr[option], SP_arr[K], SP_arr[K+1], SP_old[index], SP_old[option]);
         }
         else
         {
            //if ((index==177 && option==493) || (index==493 && option==177)){mexPrintf("2) E1K1=%e, E1K2=%e\n", E_1_K1, E_1_K2);drawnow();}
            new_ks[0] = K+1;
            new_ks[1] = K;
            max_E = move_split_calc_delta(SP_arr[index], SP_arr[option], SP_arr[K+1], SP_arr[K], SP_old[index], SP_old[option]);
         }
      }
      else
      {
         E_1_K1 = move_split_calc_delta(SP_arr[index], SP_arr[option], SP_arr[K], SP_arr[K+1], SP_old[index], SP_old[option]);
         E_1_K2 = move_split_calc_delta(SP_arr[index], SP_arr[option], SP_arr[K+1], SP_arr[K], SP_old[index], SP_old[option]);

         if (E_1_K1 > E_1_K2)
         {
            new_ks[0] = K;
            new_ks[1] = K+1;
            max_E = E_1_K1;
         }
         else
         {
            new_ks[0] = K+1;
            new_ks[1] = K;
            max_E = E_1_K2;
         }
      }
   }
}



void IMG::move_split_SP(int index)
{
   int num_SP = 2;
   if (SP_arr[index]!=NULL && !(SP_arr[index]->isempty()) && SP_arr[index]->get_N() > num_SP)
   {

      int* new_ks = new int[num_SP];
      memset(new_ks,-1,sizeof(int)*num_SP);
      int ksplit = -1;
      double max_E = -DBL_MAX;

      //only work for num_SP == 2, so far
      move_split_SP_propose(index, num_SP, -1, max_E, ksplit, new_ks);

      // update
      if(max_E>0)
      {
         //mexPrintf("split: %f,%d,%d\n",max_E,index,new_ks[0]);
         // update the labels first
         linkedListNode<int>* node = SP_arr[new_ks[0]]->pixels.getFirst();
         while (node != NULL)
         {
            label[node->getData()] = index;
            node = node->getNext();
         }

         // merge the super pixels
         SP_arr[index]->merge_with(SP_arr[new_ks[0]], label, border_ptr, xdim, ydim);
         SP_arr[new_ks[0]]->empty();
         delete SP_arr[new_ks[0]];
         SP_arr[new_ks[0]] = NULL;

         SP_changed[index] = true;
         if (ksplit<0)
         {
            // splitting into a new super pixel
            SP_changed[new_ks[1]] = true;
            // move it to the right spot
            if (new_ks[1]!=K)
            {
               node = SP_arr[new_ks[1]]->pixels.getFirst();
               while (node != NULL)
               {
                  label[node->getData()] = K;
                  node = node->getNext();
               }
               if (SP_arr[K]!=NULL)
                  mexPrintf("%d\n",SP_arr[K]->get_N());
               SP_arr[K] = SP_arr[new_ks[1]];
               SP_old[K] = SP_old[new_ks[1]];
               SP_arr[new_ks[1]] = NULL;
            }
            K++;
            max_UID++;
         }
         else
         {
            // splitting into an old super pixel
            SP_changed[ksplit] = true;
            node = SP_arr[new_ks[1]]->pixels.getFirst();
            while (node != NULL)
            {
               label[node->getData()] = ksplit;
               node = node->getNext();
            }
            SP_arr[ksplit]->merge_with(SP_arr[new_ks[1]], label, border_ptr, xdim, ydim);
            SP_arr[new_ks[1]]->empty();
            delete SP_arr[new_ks[1]];
            SP_arr[new_ks[1]] = NULL;
            alive_dead_changed = true;
         }
      }
      else if (ksplit!=-2)
      {
         //Recover previous label
         linkedListNode<int> *node = SP_arr[K]->pixels.getFirst();
         while (node!=NULL)
         {
            label[node->getData()] = index;
            node = node->getNext();
         }
         node = SP_arr[K+1]->pixels.getFirst();
         while (node!=NULL)
         {
            label[node->getData()] = index;
            node = node->getNext();
         }

         SP_changed[index] = false;

         for( int i = 0; i < num_SP; i++ )
         {
            SP_arr[index]->merge_with(SP_arr[K+i], label, border_ptr, xdim, ydim);
            SP_arr[K+i]->empty();
            delete SP_arr[K+i];
            SP_arr[K+i] = NULL;
         }
      }
      else
      {
         SP_changed[index] = false;
      }

      delete[] new_ks;
   }
}






// this is split stuff!
double IMG::U_dist(int index1 , int index2)
{
   double dist = 0;
   int ind1=index1*5,ind2=index2*5;
   //printf("dis:%d,%d\n",ind1,ind2);
   //x,y
   //arr(double) pos_Delta = new_pos.get_Delta();
   //arr(double) app_Delta = new_app.get_Delta();
   for( int n = 0; n < 2; n++ )
   {
     dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
     // dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
   }
   //L,a,b
   for( int n = 2; n < 5; n++ )
   {
      dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
      //dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
   }
      //printf("wa: %f\n",(dist/10000.0));
   return dist/10000.0;
}
double IMG::U_dist(int index1 , double* center)
{
   double dist = 0;
   int ind1=index1*5;

   //arr(double) pos_Delta=this->new_pos.get_Delta();
   //arr(double) app_Delta=this->new_app.get_Delta();

   for( int n = 0; n < 2; n++ ){
      dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
      //dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
   }
   //L,a,b
   for( int n = 2; n < 5; n++ ){
      dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
      //dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
   }
      //printf("wa2 : %f\n",(dist/10000.0));
      return dist/10000.0;
}
bool IMG::U_Kmeans_plusplus(int index, int *bbox, int num_SP, int numiter, int index2)
{
   bool broken = false;

   if (num_SP!=2)
      mexErrMsgTxt("Trying to split into more than 2");
   //void IMG::U_Kmeans_plusplus( int num_SP, int numiter){
   int num_pix = SP_arr[index]->get_N();
   //int num_pix = width*height;
   double* distvec = new double[num_pix];
   int* klabels = new int[num_pix];
   int* SP_sz = new int[num_SP];
   double * center = new double[num_SP * 5];
   double dist,tmp_dist;
   int change;

   //1. kmeans ++ initialization
   //first cener pt

   int seed;


   linkedListNode<int>* node = SP_arr[index]->pixels.getFirst();
   int count;
   int tmp_pos;
   while (node!=NULL)
   {
      //distvec[n++] = U_dist(seed,node->getData());
      //distvec[n++] = U_dist(node->getData(), center);
      // get the Bounding Box of SP_arr[index]
      // xmin,xmax,ymin,ymax
      // x
      tmp_pos = node->getData()%xdim;
      if(tmp_pos <bbox[0])
         bbox[0] = tmp_pos;
      if(tmp_pos >bbox[1])
         bbox[1] = tmp_pos;
      // y
      tmp_pos = node->getData()/xdim;
      if(tmp_pos <bbox[2])
         bbox[2] = tmp_pos;
      if(tmp_pos >bbox[3])
         bbox[3] = tmp_pos;
      node = node->getNext();
   }

   bool old_split1 = SP_old[index];
   bool old_split2 = (index2!=-1 && SP_old[index2]);
   if (!old_split1 && !old_split2) // both new
   {
      // first is new
      seed = SP_arr[index]->pixels.dataAt(rand()%num_pix);
      if(label[seed]!=index)
      {
         mexPrintf("SP_arr[index]->get_N()=%d\n", SP_arr[index]->get_N());
         mexErrMsgTxt("inconsistency about cluster label\n");
      }
      for( int n = 0; n < 5; n++ )
         center[n] = data[seed*5+n];

      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center);
         node = node->getNext();
      }

      // second is new
      seed = SP_arr[index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[5+m] = data[seed*5+m];
   }
   else if (old_split1 && !old_split2) // old new
   {
      arr(double) mean_pos = SP_arr[index]->get_mean_pos();
      arr(double) mean_app = SP_arr[index]->get_mean_app();
      center[0] = mean_pos[0];
      center[1] = mean_pos[1];
      center[2] = mean_app[0];
      center[3] = mean_app[1];
      center[4] = mean_app[2];

      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center);
         node = node->getNext();
      }

      // second is new
      seed = SP_arr[index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[5+m] = data[seed*5+m];
   }
   else if (!old_split1 && old_split2) // new old
   {
      arr(double) mean_pos = SP_arr[index2]->get_mean_pos();
      arr(double) mean_app = SP_arr[index2]->get_mean_app();
      center[5] = mean_pos[0];
      center[6] = mean_pos[1];
      center[7] = mean_app[0];
      center[8] = mean_app[1];
      center[9] = mean_app[2];

      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center+5);
         node = node->getNext();
      }

      // first is new
      seed = SP_arr[index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[m] = data[seed*5+m];
   }
   else // old old
   {
      arr(double) mean_pos = SP_arr[index]->get_mean_pos();
      arr(double) mean_app = SP_arr[index]->get_mean_app();
      center[0] = mean_pos[0];
      center[1] = mean_pos[1];
      center[2] = mean_app[0];
      center[3] = mean_app[1];
      center[4] = mean_app[2];

      mean_pos = SP_arr[index2]->get_mean_pos();
      mean_app = SP_arr[index2]->get_mean_app();
      center[5] = mean_pos[0];
      center[6] = mean_pos[1];
      center[7] = mean_app[0];
      center[8] = mean_app[1];
      center[9] = mean_app[2];
   }


   //2. kmeans ++ iterations
   change = 0;
   for( int itr = 0; itr < numiter; itr++ ){
      memset(distvec, 0x7f, sizeof(double)*num_pix);
      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {

         for( int i = 0; i < num_SP; i++ ){
            tmp_dist = U_dist(node->getData(),&center[i*5]);
            //				printf("%d,%d,%f\n",n,node->getData(),tmp_dist);
            if(tmp_dist < distvec[count] ){
               distvec[count] = tmp_dist;
               klabels[count]  = i;
               change = 1;
            }
         }
         node = node->getNext();
         count++;
      }

      if (change ==0){
         //no change happened... Kmeans totally stuck
         break;
      }
      memset(center,0,sizeof(double)*num_SP*5);
      memset(SP_sz,0,sizeof(int)*num_SP);


      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         //printf("dist %d,%d,%d,%d\n",r,c,counter,klabels[counter]);
         // klabels[n]==0 && old_split1 then don't update
         // klabels[n]==1 && old_split2 then don't update
         // if (klabels[n]!=old_split1-1 && klabels[n]!=old_split2-1)
         if ((klabels[count]!=0 || !old_split1) && (klabels[count]!=1 || !old_split2))
            for( int j = 0; j < 5; j++ )
               center[klabels[count]*5+j] += data[(node->getData())*5+j];
         SP_sz[klabels[count]] += 1;
         count++;
         node = node->getNext();
      }

      double inv;
      for( int k = 0; k < num_SP; k++ )
      {
         if (SP_sz[k]>0)
         {
            if ((k==0 && !old_split1) || (k==1 && !old_split2))
            {
               inv = 1.0/double(SP_sz[k]);
               for( int kk = 0; kk < 5; kk++ )
                  center[k*5+kk] *= inv;
            }
         }
         else
         {
            if (old_split1 || old_split2)
               broken = true;
            else
               mexErrMsgTxt("one cluster removed... shouldn't happen\n");
         }
      }
   }

   if (!broken)
   {
      // change label accordingly
      node = SP_arr[index]->pixels.getFirst();
      count = 0;
      while (node != NULL)
      {
         label[node->getData()] = K+klabels[count++];
         node = node->getNext();
      }
   }

   delete[] distvec;
   delete[] SP_sz;
   delete[] center;
   delete[] klabels;
   return broken;
}



// KMerge Version 2: grow region

void IMG::U_connect_newSP(int *bbox,int num_SP){
   int b, c;
   int min_x = bbox[0];
   int min_y = bbox[2];
   int max_x = bbox[1];
   int max_y = bbox[3];
   int tmp_xdim = max_x-min_x+1;
   int tmp_ydim = max_y-min_y+1;
   int tmp_N = tmp_ydim*tmp_xdim;
   //printf("%d,%d,%d,%d\n",bbox[0],bbox[1],bbox[2],bbox[3]);
   const int dx4[4] = {-1,  0,  1,  0};
   const int dy4[4] = { 0, -1,  0,  1};

   arr(int) nlabels = allocate_memory<int>(tmp_N,-1);
   int oindex = 0, adjlabel = 0, label_count =0;

   int base_ind = min_y*xdim+min_x,tmp_ind,tmp_x,tmp_y,new_ind,new_nind,new_x,new_y;
   linkedListNode<int>* node;
   linkedList<int> check_labels;
   for(int j =0;j<num_SP;j++)
      check_labels.addNodeEnd(K+j);

   for( int j = 0; j < tmp_ydim; j++ )//tmp_ydim
   {
      for( int k = 0; k < tmp_xdim; k++ )//tmp_xdim
      {

         tmp_ind = base_ind+k;
         if(U_FindArray(label[tmp_ind],check_labels)!=-1 && nlabels[oindex]==-1)
         {
         //printf("do it %d,%d,%d,%d,%d,%d\n",j,k,tmp_ind,base_ind,min_y*xdim,min_x);
            // first find an unlabeled pixel in labels
            nlabels[oindex] = label_count;
            //--------------------
            // Start a new segment of index oindex
            //--------------------
            if(SP_arr[K + label_count]!= NULL)
               mexErrMsgTxt("SP should be null...\n");
            //mexErrMsgTxt(strcat(,"th SP should be null...\n"));
            SP_arr[K + label_count] = new SP(new_pos, new_app, max_UID);
            //can't decide whether it will be border yet
            SP_arr[K + label_count]->add_pixel(data, tmp_ind, false, pixel_ptr[tmp_ind], border_ptr[tmp_ind], boundary_mask[tmp_ind]);

            node = SP_arr[K + label_count]->pixels.getFirst();
            //int cc=1;
            while (node!=NULL)
            {
               tmp_x = node->getData()%xdim;
               tmp_y = node->getData()/xdim;
               for( int n = 0; n < 4; n++ )
               {
                  new_x = tmp_x+dx4[n];
                  new_y = tmp_y+dy4[n];
                  if( (new_x >= min_x && new_x <= max_x) && (new_y >= min_y && new_y <= max_y) )
                  {
                     new_ind = new_y*xdim + new_x;
                     new_nind = (new_x-min_x)+(new_y-min_y)*tmp_xdim;
                     if(nlabels[new_nind] ==-1 && label[new_ind] == label[tmp_ind] )
                     {
                        // should always update labels before adding pixel, otherwise
                        // U_check_border_pix will be wrong!
                        nlabels[new_nind] = label_count;
                        label[new_ind] = K+label_count;
                        SP_arr[K + label_count]->add_pixel(data, new_ind, false, pixel_ptr[new_ind], border_ptr[new_ind], boundary_mask[new_ind]);
                     }
                  }
               }
               node = node->getNext();
            }
            label[tmp_ind] = K+label_count;

            label_count++;
         }
         oindex ++;
      }
      base_ind += xdim;
   }

   // Now is the time to clean up the border pixels and neighbor ids
   for (int k=K; k<K+label_count; k++){
      SP_arr[k]->fix_borders(label, border_ptr, xdim, ydim);
      U_fix_neighbors_self(k);
   }

   if (label_count>num_SP)
   {
      //printf("oops... kmeans gives %d connected components\n",label_count);
      int* pix_counts =  new int[label_count];
      int* ind_counts =  new int[label_count];
      for(int i =0; i<label_count; i++){
         pix_counts[i] = SP_arr[K+i]->N;
         ind_counts[i] = i;
            if(i>=num_SP) check_labels.addNode(K+i);
            //printf("SP: %d, %d pixels with label %d \n",K+i,SP_arr[K+i]->N,label[SP_arr[K+i]->pixels.getFirst()->getData()]);
      }
      U_quicksort<int>(pix_counts,ind_counts,0,label_count-1);

      //need to merge smallest one to its best neighbour
         //if(K==1719)
      //	mexErrMsgTxt("the redundant SP doesn't merge within the range...\n");

      double max_E;
      int max_k;
      arr(bool) neighbors = allocate_memory<bool>(K+label_count,false);
      for(int i =0; i<label_count-num_SP; i++)
      {
         max_k = -1;
         max_E = -DBL_MAX;
         move_merge_SP_propose_region( K+ind_counts[i], neighbors,check_labels,max_E,max_k);
         //printf("to merge %d,%d,%d\n",i,K+ind_counts[i],max_k);

         if(max_k==-1||max_E==-1)
            mexErrMsgTxt("the redundant SP doesn't merge within the range...\n");
         node = SP_arr[K+ind_counts[i]]->pixels.getFirst();
         while(node !=NULL){
                  //if(K==689) printf("mmmerme:%d,%d\n",label[node->getData()],label[SP_arr[max_k]->pixels.getFirst()->getData()]);
            label[node->getData()] = max_k;
            node = node->getNext();
         }
         SP_arr[max_k]->merge_with(SP_arr[K+ind_counts[i]],label, border_ptr, xdim, ydim);
         delete SP_arr[K+ind_counts[i]];
         SP_arr[K+ind_counts[i]] = NULL;
      }

      // relabel SPs
      //printf("relabel\n");
      max_k = K;
      for(int i =0; i<num_SP; i++)
      {
         while(SP_arr[max_k]==NULL)
         {
            max_k++;
            if(max_k>N)
               mexErrMsgTxt("should be more nonempty SPs...\n");
         }
         if(max_k!=K+i)
         {
            //fetch stuff far behind to here
            SP_arr[K+i] = SP_arr[max_k];
            SP_arr[max_k] = NULL;
            //printf("pair up %d,%d\n",(K+i),max_k);
            //relabel to the first two ...
            node = SP_arr[K+i]->pixels.getFirst();
            while(node !=NULL)
            {
               label[node->getData()] = K+i;
               node = node->getNext();
            }
         }
         max_k++;
      }

      delete[] pix_counts;
      delete[] ind_counts;
      deallocate_memory(neighbors);
   }
   deallocate_memory(nlabels);
}



