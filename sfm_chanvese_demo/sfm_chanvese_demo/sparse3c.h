#ifndef __SPARSE3C_H
#define __SPARSE3C_H

#include <stdlib.h>
#include "mex.h"
#include <matrix.h>
#include <math.h>
#include "llist.h"
#include "lsops3c.h"
#include "energy3c.h"

#define DIMX dims[1]
#define DIMY dims[0]
#define DIMZ dims[2]
#define DIMXY dims[3]
#define NUMEL dims[4]

#define OFFX dims[0]
#define OFFY 1
#define OFFZ dims[3]

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void show_countdown(int iter,int i, int *countdown, mxArray **ppphi);

mxArray *prep_C_output(LL *Lz,long *dims, double *phi);

// "Soft Plaque Detection and Vessel Segmentation" 
// - Lankton, et. al.
void lrbac_vessel_yz(double *img, double *phi, double *label, long *dims,
                  LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                  int iter, double rad, double lambda, double dthresh, mxArray **plhs,
                  int display);

// "Soft Plaque Detection and Vessel Segmentation" 
// - Lankton, et. al.
void lrbac_vessel_cv(double *img, double *phi, double *label, long *dims,
                  LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                  int iter, double rad, double lambda, double dthresh, mxArray **plhs,
                  int display);

// "Localizing Region-based Active Contours" 
// - Lankton, et. al.
void lrbac_chanvese(double *img, double *phi, double *label, long *dims,
                    LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                    int iter, double rad, double lambda,  mxArray **plhs,int display);

// "Image Segmentation Using Active Contours Driven by the Bhattacharyya
// Gradient Flow" - Michailovich, et al.
void bhattacharyya(double *img, double *phi, double *label, long *dims,
                   LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
                   int iter,double lambda,  mxArray **plhs,int display);

// "Active Contours Without Edges" - Chan and Vese
void chanvese(double *img, double *phi, double *label, long *dims,
              LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
              int iter,double lambda,  mxArray **plhs,int display);

// "A Variational Framework for Active and Adaptive Segmentation of Vector
// Valued Images" - Rousson and Deriche
void meanvar(double *img, double *phi, double *label, long *dims,
             LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
             int iter,double lambda,  mxArray **plhs,int display);

// "A Fully Global Approach to Image Segmentation via Coupled Curve Evolution
// Equations" - Yezzi, et. al.
void yezzi(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,  mxArray **plhs,int display);

void grow(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,  mxArray **plhs,int display);

// "Soft Plaque Detection and Vessel Segmentation:
// - Lankton et. al.
void shrink(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
            int iter,double rad,double lambda, mxArray **plhs,int display);

void kappa(double *img, double *phi, double *label, long *dims,
           LL *Lz, LL *Ln1, LL *Lp1, LL *Ln2, LL *Lp2, LL *Lin2out, LL *Lout2in,
           int iter,double lambda,  mxArray **plhs,int display);

#endif
