/*********************************************************************
 * sfm_local_chanvese_mex.cpp
 *
 * This file performs the statistical computations needed for running 
 * the localized chan-vese active contour segmentation energy using the 
 * Sparse Field method presented by Whitaker.
 * 
 * written by: Shawn Lankton (4/17/2009) - www.shawnlankton.com
 *
 ********************************************************************/

#include "sparse3c.h"

//[phi C L] = ls_sparse(img,mask,iterations,display,lam,rad,dthresh);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  //declare variables
  mxArray *phi_out, *C_out, *label_out, *mxPhi, *C_in;
  const mxArray *mxImg;

  const mwSize *mdims;
  double *img, *phi, *B, *mask, *C, *label;
  double *F;
  double usum, vsum, dthresh, lambda, rad;
  int    iter,countdown,display;
  long    dims[5];
  long    dimx, dimy, dimz, numdims;
  LL *Lz, *Ln1, *Ln2, *Lp1, *Lp2;
  LL *Sz, *Sn1, *Sn2, *Sp1, *Sp2;
  LL *Lin2out, *Lout2in;

  //figure out dimensions
  mdims = mxGetDimensions(prhs[0]);
  dims[2] = 1; dims[1] = 1;
  numdims = mxGetNumberOfDimensions(prhs[0]);
  switch(numdims){
  case 3: dimz = (int)mdims[2]; dims[2] = dimz;
  case 2: dimx = (int)mdims[1]; dims[1] = dimx;
  case 1: dimy = (int)mdims[0]; dims[0] = dimy;}
  dims[3] = dims[0]*dims[1]; dims[4] = dims[0]*dims[1]*dims[2];


  // associate inputs;
  img     = mxGetPr(prhs[0]);
  mask    = mxGetPr(prhs[1]);
  iter    = (int)mxGetPr(prhs[2])[0];
  if(nrhs>3) display= (int)(mxGetPr(prhs[3])[0]);     else display = 1;
  if(nrhs>4) lambda = (double)(mxGetPr(prhs[4])[0]);  else lambda = .1;
  if(nrhs>5) rad    = (double)(mxGetPr(prhs[5])[0]);  else rad = 10;
  if(nrhs>6) dthresh= (double)(mxGetPr(prhs[6])[0]);  else dthresh = 500;

  //associate outputs;
  phi_out  = plhs[0] = mxDuplicateArray(prhs[1]);
  label_out= plhs[2] = mxDuplicateArray(prhs[1]);
  mxPhi = plhs[0];
  phi  = mxGetPr(phi_out);
  // C is not setup until after the segmentation
  label= mxGetPr(label_out);


  //create linked lists
  Lz  = ll_create();
  Ln1 = ll_create();
  Ln2 = ll_create();
  Lp1 = ll_create();
  Lp2 = ll_create();
  Lin2out = ll_create();
  Lout2in = ll_create();

  //initialize lists, phi, and labels
  ls_mask2phi3c(mask,phi,label,dims,Lz,Ln1,Ln2,Lp1,Lp2);

  //------------------------------
  // IMPLEMENTED ENERGIES
  //-----------------------------
  //lrbac_vessel_yz(img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,rad,lambda,dthresh,plhs,display);
  //lrbac_vessel_cv(img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,rad,lambda,dthresh,plhs,display);
  lrbac_chanvese( img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,rad,lambda,plhs,display);
  //bhattacharyya(  img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);
  //chanvese(       img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);
  //meanvar(        img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);
  //yezzi(          img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);
  //grow(           img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);
  //shrink(         img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,rad,lambda,plhs,display);
  //kappa(          img,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in,iter,lambda,plhs,display);


  //prepare "list" output
  plhs[1] = prep_C_output(Lz,dims,phi);

  //destroy linked lists
  ll_destroy(Lz);
  ll_destroy(Ln1);
  ll_destroy(Ln2);
  ll_destroy(Lp1);
  ll_destroy(Lp2);
  ll_destroy(Lin2out);
  ll_destroy(Lout2in);
  return;
}


void lrbac_vessel_yz(double *img, double *phi, double *label, long *dims,
                     LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                     int iter, double rad, double lambda, double dthresh, mxArray **plhs,
                     int display){
  double *F;
  double scale[1]; scale[0]=0;
  int countdown;

  //initialize datastructures and statistics
  en_lrbac_init(Lz,img,phi,dims,rad);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_lrbac_vessel_yz_compute(Lz,phi,img,dims, scale,lambda,rad,dthresh);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_lrbac_vessel_yz_update(img, dims, Lin2out, Lout2in, rad,dthresh);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
  //destroy old datastructures
  en_lrbac_destroy();
}

void lrbac_vessel_cv(double *img, double *phi, double *label, long *dims,
                     LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                     int iter, double rad, double lambda, double dthresh, mxArray **plhs,
                     int display){
  double *F;
  double scale[1]; scale[0]=0;
  int countdown;

  //initialize datastructures and statistics
  en_lrbac_init(Lz,img,phi,dims,rad);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_lrbac_vessel_cv_compute(Lz,phi,img,dims, scale,lambda,rad,dthresh);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_lrbac_vessel_cv_update(img, dims, Lin2out, Lout2in, rad,dthresh);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
  //destroy old datastructures
  en_lrbac_destroy();
}

void lrbac_chanvese(double *img, double *phi, double *label, long *dims,
                    LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                    int iter, double rad, double lambda,  mxArray **plhs,int display){
  double *F;
  double scale[1]; scale[0]=0;
  int countdown;

  //initialize datastructures and statistics
  en_lrbac_init(Lz,img,phi,dims,rad);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_lrbac_compute(Lz,phi,img,dims, scale,lambda,rad);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_lrbac_update(img, dims, Lin2out, Lout2in, rad);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
  //destroy old datastructures
  en_lrbac_destroy();
}

void chanvese(double *img, double *phi, double *label, long *dims,
              LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
              int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  double scale[1]; scale[0] = 0;
  int countdown;

  //initialize datastructures and statistics
  en_chanvese_init(img,phi,dims);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_chanvese_compute(Lz,phi,img,dims,scale,lambda);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_chanvese_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
}

void meanvar(double *img, double *phi, double *label, long *dims,
             LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
             int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  double scale[1]; scale[0] = 0;
  int countdown;

  //initialize datastructures and statistics
  en_meanvar_init(img,phi,dims);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_meanvar_compute(Lz,phi,img,dims,scale,lambda);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_meanvar_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
}

void bhattacharyya(double *img, double *phi, double *label, long *dims,
                   LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                   int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  double scale[1]; scale[0] = 0;
  int countdown;

  //initialize datastructures and statistics
  en_bhattacharyya_init(img,phi,dims);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_bhattacharyya_compute(Lz,phi,img,dims,scale,lambda);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_bhattacharyya_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
  en_bhattacharyya_destroy();
}


void yezzi(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  double scale[1]; scale[0]=0;
  int countdown;

  //initialize datastructures and statistics
  en_yezzi_init(Lz,img,phi,dims);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_yezzi_compute(Lz,phi,img,dims,scale,lambda);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_yezzi_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
}

void grow(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
          int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  int countdown;

  for(int i=0;i<iter;i++){
    //compute force
    F = en_grow_compute(Lz, img, phi,dims,lambda, 800);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_null_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
}

void shrink(double *img, double *phi, double *label, long *dims,
            LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
            int iter,double rad,double lambda, mxArray **plhs,int display){
  double *F;
  int countdown;
  double scale[1]; scale[0] = 0;

  en_lrbac_init(Lz,img,phi,dims,rad);
  for(int i=0;i<iter;i++){
    //compute force
    F = en_shrink_compute(Lz,img,phi,dims,rad,lambda,scale);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_lrbac_update(img, dims, Lin2out, Lout2in, rad);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
  //destroy old datastructures
  en_lrbac_destroy();
}

void kappa(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,  mxArray **plhs,int display){
  double *F;
  int countdown;

  for(int i=0;i<iter;i++){
    //compute force
    F = en_kappa_compute(Lz,phi,dims);
    //perform iteration
    ls_iteration(F,phi,label,dims,Lz,Ln1,Lp1,Ln2,Lp2,Lin2out,Lout2in);
    //update statistics
    en_null_update(img, dims, Lin2out, Lout2in);

    //display stuff (maybe)
    if(display) show_countdown(iter,i,&countdown,plhs);
  }
}

mxArray *prep_C_output(LL *Lz,long *dims,double *phi){
  mxArray *out;
  double *data, *C;
  int n = 0;
  int x,y,z,idx;

  if(Lz == NULL) return mxCreateDoubleMatrix(1,1,mxREAL);
  out = mxCreateDoubleMatrix(Lz->length,1,mxREAL);
  C = mxGetPr(out);

  ll_init(Lz);
  while(Lz->curr != NULL){
    C[n] =(double)(Lz->curr->x*DIMY+Lz->curr->y+Lz->curr->z*DIMXY);
    ll_step(Lz); n++;
  }
  return out;
}

void show_countdown(int iter,int i,int *countdown, mxArray **ppphi){
  double loopmax;
  int period = 10;
  if(countdown[0] <= 0){
    countdown[0] = iter/period;
    loopmax = (int)((double)i/(double)iter*(double)period)*(100/period);
    if(i>0){
      mexCallMATLAB(0,NULL,1,ppphi,"visualize_phi");
      mexPrintf(" %d%% complete\n",(int)loopmax);
    }
  }
  else{
    countdown[0]--;
  }
  if(i==iter-1){
    mexCallMATLAB(0,NULL,1,ppphi,"visualize_phi");
    mexPrintf("100%% complete\n");
  }
}
