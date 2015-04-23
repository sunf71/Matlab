#include "lsops3c.h"

void ls_iteration(double *F, double *phi, double* label, long* dims, 
                  LL* Lz, LL* Ln1, LL* Lp1, LL *Ln2, LL *Lp2, 
                  LL *Lin2out, LL* Lout2in){
  int x,y,z,i,idx;
  int u,d,r,l,f,b;
  double p, phi_old;
  LL *Sz, *Sn1, *Sp1, *Sn2, *Sp2;

  // create 'changing status' lists
  Sz  = ll_create();
  Sn1 = ll_create();
  Sp1 = ll_create();
  Sn2 = ll_create();
  Sp2 = ll_create();

  // #1) Normalize F
  double Fmax = .001;
  for(i=0;i<Lz->length;i++){
    if(fabs(F[i])>Fmax) Fmax = fabs(F[i]);
  }

  for(i=0;i<Lz->length;i++){ 
    F[i] = F[i]/Fmax*0.4;
  }

  // #2) add F to phi(Lz), create Sn1 & Sp1 
  //                                             ========
  //     (a) scan Lz values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lz); i = 0;
  while(Lz->curr != NULL){
    x= Lz->curr->x; y= Lz->curr->y; z= Lz->curr->z; idx= Lz->curr->idx;
    phi_old = phi[idx];

    phi[idx] = phi[idx]+F[i];

    //check to see if point crossed interface
    if(phi_old<=0 && phi[idx]>0 ){
      ll_pushnew(Lin2out,x,y,z,idx);
    }
    if(phi_old>0  && phi[idx]<=0){
      ll_pushnew(Lout2in,x,y,z,idx);
    }
    
    if(phi[idx] > .5){
      ll_push(Sp1, ll_remcurr(Lz)); 
    }
    else if(phi[idx] < -.5){
      ll_push(Sn1, ll_remcurr(Lz));
    }
    else{
      ll_step(Lz);
    }
    i++; //increment index into F
  }
  if(F!= NULL) mxFree(F); // free F (no longer needed);

  // #3) update Ln1,Ln2,Lp1,Lp2
  //                                    ==========
  //     (c) scan Ln1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Ln1);
  while(Ln1->curr != NULL){
    x= Ln1->curr->x; y= Ln1->curr->y; z= Ln1->curr->z; idx= Ln1->curr->idx;
    p = ls_max_hood_onlevel(idx, x, y, z, dims, phi,label,0);
    if(p>=-0.5){        // found something
      phi[idx] = p-1;

      if(phi[idx]>=-0.5){
        ll_push(Sz,ll_remcurr(Ln1));
      }
      else if(phi[idx]<-1.5){
        ll_push(Sn2,ll_remcurr(Ln1));
      }
      else ll_step(Ln1);
    }
    else{
      ll_push(Sn2,ll_remcurr(Ln1));
    }
  }

  //                                                      ========
  //     (c) scan Lp1 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lp1);
  while(Lp1->curr != NULL){
    x= Lp1->curr->x; y= Lp1->curr->y; z= Lp1->curr->z; idx= Lp1->curr->idx;
    p = ls_min_hood_onlevel(idx, x, y, z, dims, phi,label,0);
    if(p<=0.5){         // found something
      phi_old = phi[idx];
      phi[idx] = p+1;

      if(phi[idx]<=0.5){
        ll_push(Sz,ll_remcurr(Lp1));
      }
      else if(phi[idx]>1.5){
        ll_push(Sp2,ll_remcurr(Lp1));
      }
      else ll_step(Lp1);
    }
    else{
      ll_push(Sp2,ll_remcurr(Lp1));
    }
  }

  //                         ===========
  //     (c) scan Ln2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Ln2);
  while(Ln2->curr != NULL){
    x = Ln2->curr->x; y = Ln2->curr->y; z = Ln2->curr->z; idx= Ln2->curr->idx;
    p = ls_max_hood_onlevel(idx, x, y, z, dims, phi,label,-1);
    if(p>=-1.5){         // found something
      phi[idx] = p-1;
      if(phi[idx]>=-1.5){
        ll_push(Sn1,ll_remcurr(Ln2));
      }
      else if(phi[idx]<-2.5){
        ll_remcurr_free(Ln2);
        phi[idx] = -3; label[idx] = -3;
      }
      else ll_step(Ln2);
    }
    else{
      ll_remcurr_free(Ln2);
      phi[idx] = -3; label[idx] = -3;
    }
  }

  //                                                              =========
  //     (d) scan Lp2 values [-2.5 -1.5)[-1.5 -.5)[-.5 .5](.5 1.5](1.5 2.5]
  ll_init(Lp2);
  while(Lp2->curr != NULL){
    x = Lp2->curr->x; y = Lp2->curr->y; z = Lp2->curr->z; idx= Lp2->curr->idx;
    p = ls_min_hood_onlevel(idx, x, y, z, dims, phi,label,1);
    if(p<=1.5){         // found something
      phi[idx] = p+1;
      if(phi[idx]<=1.5){
        ll_push(Sp1,ll_remcurr(Lp2));
      }
      else if(phi[idx]>2.5){
        ll_remcurr_free(Lp2);
        phi[idx] = 3; label[idx] = 3;
      }
      else ll_step(Lp2);
    }
    else{
      ll_remcurr_free(Lp2);
      phi[idx] = 3; label[idx] = 3;
    }
  }

  // #4) Deal with S-lists Sz,Sn1,Sp1,Sn2,Sp2
  //     (a) Scan Sz
  ll_init(Sz);
  while(Sz->curr != NULL){
    idx= Sz->curr->idx;
    ll_push(Lz,ll_remcurr(Sz));
    label[idx] = 0;
  }

  //     (b) Scan Sn1
  ll_init(Sn1);
  while(Sn1->curr != NULL){
    x = Sn1->curr->x; y = Sn1->curr->y; z = Sn1->curr->z; idx = Sn1->curr->idx;
    ll_push(Ln1,ll_remcurr(Sn1));
    label[idx] = -1;
    if(((y+1)<DIMY) && phi[idx+OFFY]==-3){
      ll_pushnew(Sn2,x,y+1,z,idx+OFFY); 
      phi[idx+OFFY] = phi[idx] - 1; 
    }//up
    if(((y-1)>=0)   && phi[idx-OFFY]==-3){
      ll_pushnew(Sn2,x,y-1,z,idx-OFFY);
      phi[idx-OFFY] = phi[idx] - 1; 
    }//down
    if(((x+1)<DIMX) && phi[idx+OFFX]==-3){
      ll_pushnew(Sn2,x+1,y,z,idx+OFFX);
      phi[idx+OFFX] = phi[idx] - 1; 
    }//right
    if(((x-1)>=0)   && phi[idx-OFFX]==-3){
      ll_pushnew(Sn2,x-1,y,z,idx-OFFX);
      phi[idx-OFFX] = phi[idx] - 1; 
    }//left
    if(((z+1)<DIMZ) && phi[idx+OFFZ]==-3){
      ll_pushnew(Sn2,x,y,z+1,idx+OFFZ);
      phi[idx+OFFZ] = phi[idx] - 1; 
    }//front
    if(((z-1)>=0)   && phi[idx-OFFZ]==-3){
      ll_pushnew(Sn2,x,y,z-1,idx-OFFZ);
      phi[idx-OFFZ] = phi[idx] - 1; 
    }//back
  }

  //     (c) Scan Sp1
  ll_init(Sp1);
  while(Sp1->curr != NULL){
    x = Sp1->curr->x; y = Sp1->curr->y; z = Sp1->curr->z; idx=Sp1->curr->idx;
    ll_push(Lp1,ll_remcurr(Sp1));
    label[idx] = 1;
    if(((y+1)<DIMY) && phi[idx+OFFY]==3){
      ll_pushnew(Sp2,x,y+1,z,idx+OFFY); 
      phi[idx+OFFY] = phi[idx] + 1;
    }//up
    if(((y-1)>=0)   && phi[idx-OFFY]==3){
      ll_pushnew(Sp2,x,y-1,z,idx-OFFY); 
      phi[idx-OFFY] = phi[idx] + 1;
    }//down
    if(((x+1)<DIMX) && phi[idx+OFFX]==3){
      ll_pushnew(Sp2,x+1,y,z,idx+OFFX); 
      phi[idx+OFFX] = phi[idx] + 1;
    }//right
    if(((x-1)>=0)   && phi[idx-OFFX]==3){
      ll_pushnew(Sp2,x-1,y,z,idx-OFFX); 
      phi[idx-OFFX] = phi[idx] + 1;
    }//left
    if(((z+1)<DIMZ) && phi[idx+OFFZ]==3){
      ll_pushnew(Sp2,x,y,z+1,idx+OFFZ); 
      phi[idx+OFFZ] = phi[idx] + 1;
    }//front
    if(((z-1)>=0)   && phi[idx-OFFZ]==3){
      ll_pushnew(Sp2,x,y,z-1,idx-OFFZ); 
      phi[idx-OFFZ] = phi[idx] + 1;
    }//back
  }

  //     (d) Scan Sn2
  ll_init(Sn2);
  while(Sn2->curr != NULL){
    idx = Sn2->curr->idx;
    ll_push(Ln2,ll_remcurr(Sn2));
    label[idx] = -2;
  }

  //     (e) Scan Sp2
  ll_init(Sp2);
  while(Sp2->curr != NULL){
    idx = Sp2->curr->idx;
    ll_push(Lp2,ll_remcurr(Sp2));
    label[idx] = 2;
  }

  ll_destroy(Sz);
  ll_destroy(Sn1);
  ll_destroy(Sp1);
  ll_destroy(Sn2);
  ll_destroy(Sp2);
}

double ls_max_hood_onlevel(int idx, long x, long y, long z, long *dims, double *phi, double *label, double level){
  double pmax = -3;
  if(((y+1)<DIMY) && (label[idx+OFFY]>=level) && (phi[idx+OFFY]>pmax)) pmax = phi[idx+OFFY]; 
  if(((y-1)>=0  ) && (label[idx-OFFY]>=level) && (phi[idx-OFFY]>pmax)) pmax = phi[idx-OFFY]; 
  if(((x+1)<DIMX) && (label[idx+OFFX]>=level) && (phi[idx+OFFX]>pmax)) pmax = phi[idx+OFFX]; 
  if(((x-1)>=0  ) && (label[idx-OFFX]>=level) && (phi[idx-OFFX]>pmax)) pmax = phi[idx-OFFX]; 
  if(((z+1)<DIMZ) && (label[idx+OFFZ]>=level) && (phi[idx+OFFZ]>pmax)) pmax = phi[idx+OFFZ]; 
  if(((z-1)>=0  ) && (label[idx-OFFZ]>=level) && (phi[idx-OFFZ]>pmax)) pmax = phi[idx-OFFZ]; 
  return pmax;
}

double ls_min_hood_onlevel(int idx, long x, long y, long z, long *dims, double *phi, double *label, double level){
  double pmin = 3;
  if(((y+1)<DIMY) && (label[idx+OFFY]<=level) && (phi[idx+OFFY]<pmin)) pmin = phi[idx+OFFY]; 
  if(((y-1)>=0  ) && (label[idx-OFFY]<=level) && (phi[idx-OFFY]<pmin)) pmin = phi[idx-OFFY]; 
  if(((x+1)<DIMX) && (label[idx+OFFX]<=level) && (phi[idx+OFFX]<pmin)) pmin = phi[idx+OFFX]; 
  if(((x-1)>=0  ) && (label[idx-OFFX]<=level) && (phi[idx-OFFX]<pmin)) pmin = phi[idx-OFFX]; 
  if(((z+1)<DIMZ) && (label[idx+OFFZ]<=level) && (phi[idx+OFFZ]<pmin)) pmin = phi[idx+OFFZ]; 
  if(((z-1)>=0  ) && (label[idx-OFFZ]<=level) && (phi[idx-OFFZ]<pmin)) pmin = phi[idx-OFFZ]; 
  return pmin;
}


void ls_mask2phi3c(double* mask, double* phi, double* label, long* dims, LL *Lz, LL *Ln1, LL *Ln2, LL *Lp1, LL *Lp2){
  int x,y,z,idx;
  int i,j,k;
  int u,d,r,l,f,b;
  int  flag=0;

  //find 'interface' and mark as 0, create Lz
  for(x=0;x<DIMX;x++) for(y=0;y<DIMY;y++) for(z=0;z<DIMZ;z++){
    idx = (int)(z*DIMXY+x*DIMY+y); 

    //mark the inside and outside of label and phi
    if(mask[idx]==0){ label[idx] =  3; phi[idx]= 3; }
    else            { label[idx] = -3; phi[idx]=-3; }

    if(mask[idx] == 1){
      flag = 0;
      //if any neighbors are 1;
      if(((y+1)<DIMY) && mask[idx+OFFY]==0){flag = 1;}//up
      if(((y-1)>=0)   && mask[idx-OFFY]==0){flag = 1;}//down
      if(((x+1)<DIMX) && mask[idx+OFFX]==0){flag = 1;}//right
      if(((x-1)>=0)   && mask[idx-OFFX]==0){flag = 1;}//left
      if(((z+1)<DIMZ) && mask[idx+OFFZ]==0){flag = 1;}//front
      if(((z-1)>=0)   && mask[idx-OFFZ]==0){flag = 1;}//back
      if(flag){
        ll_pushnew(Lz,x,y,z,idx);
        label[idx] = 0;
        phi[idx] = 0;
      }
    }
  }

  //scan Lz to create Ln1 and Lp1
  ll_init(Lz);
  while(Lz->curr != NULL){
    x = Lz->curr->x; y = Lz->curr->y; z = Lz->curr->z; idx = Lz->curr->idx;
    
    if(((y+1)<DIMY) && abs((int)label[idx+OFFY])==3){//up
      if(phi[idx+OFFY]<0){//in
        label[idx+OFFY]=-1; phi[idx+OFFY]=-1;
        ll_pushnew(Ln1,x,y+1,z,idx+OFFY);
      }
      else{                //out
        label[idx+OFFY]=1; phi[idx+OFFY]=1;
        ll_pushnew(Lp1,x,y+1,z,idx+OFFY);
      }
    }
    if(((y-1)>=0)   && abs((int)label[idx-OFFY])==3){//down
      if(phi[idx-OFFY]<0){ //in
        label[idx-OFFY]=-1; phi[idx-OFFY]=-1;
        ll_pushnew(Ln1,x,y-1,z,idx-OFFY);
      }
      else{                //out
        label[idx-OFFY]=1; phi[idx-OFFY]=1;
        ll_pushnew(Lp1,x,y-1,z,idx-OFFY);
      }
    }
    if(((x+1)<DIMX) && abs((int)label[idx+OFFX])==3){//right
      if(phi[idx+OFFX]<0){//in
        label[idx+OFFX]=-1; phi[idx+OFFX]=-1;
        ll_pushnew(Ln1,x+1,y,z,idx+OFFX);
      }
      else{                //out
        label[idx+OFFX]=1; phi[idx+OFFX]=1;
        ll_pushnew(Lp1,x+1,y,z,idx+OFFX);
      }
    }
    if(((x-1)>=0)   && abs((int)label[idx-OFFX])==3){//left
      if(phi[idx-OFFX]<0){ //in
        label[idx-OFFX]=-1; phi[idx-OFFX]=-1;
        ll_pushnew(Ln1,x-1,y,z,idx-OFFX);
      }
      else{                //out
        label[idx-OFFX]=1; phi[idx-OFFX]=1;
        ll_pushnew(Lp1,x-1,y,z,idx-OFFX);
      }
    }
    if(((z+1)<DIMZ) && abs((int)label[idx+OFFZ])==3){//front
      if(phi[idx+OFFZ]<0){//in
        label[idx+OFFZ]=-1; phi[idx+OFFZ]=-1;
        ll_pushnew(Ln1,x,y,z+1,idx+OFFZ);
      }
      else{                //out
        label[idx+OFFZ]=1; phi[idx+OFFZ]=1;
        ll_pushnew(Lp1,x,y,z+1,idx+OFFZ);
      }
    }
    if(((z-1)>=0) && abs((int)label[idx-OFFZ])==3 ){//back
      if(phi[idx-OFFZ]<0){ //in
        label[idx-OFFZ]=-1; phi[idx-OFFZ]=-1;
        ll_pushnew(Ln1,x,y,z-1,idx-OFFZ);
      }
      else{                //out
        label[idx-OFFZ]=1; phi[idx-OFFZ]=1;
        ll_pushnew(Lp1,x,y,z-1,idx-OFFZ);
      }
    }

    ll_step(Lz);
  }


  //scan Ln1 to create Ln2
  ll_init(Ln1);
  while(Ln1->curr != NULL){
    x = Ln1->curr->x; y = Ln1->curr->y; z = Ln1->curr->z; idx = Ln1->curr->idx;
    
    if(((y+1)<DIMY) && label[idx+OFFY]==-3){//up
        label[idx+OFFY]=-2; phi[idx+OFFY]=-2;
        ll_pushnew(Ln2,x,y+1,z,idx+OFFY);
    }
    if(((y-1)>=0) && label[idx-OFFY]==-3){//down
        label[idx-OFFY]=-2; phi[idx-OFFY]=-2;
        ll_pushnew(Ln2,x,y-1,z,idx-OFFY);
    }
    if(((x+1)<DIMX) && label[idx+OFFX]==-3){//right
        label[idx+OFFX]=-2; phi[idx+OFFX]=-2;
        ll_pushnew(Ln2,x+1,y,z,idx+OFFX);
    }
    if(((x-1)>=0) && label[idx-OFFX]==-3){//left
        label[idx-OFFX]=-2; phi[idx-OFFX]=-2;
        ll_pushnew(Ln2,x-1,y,z,idx-OFFX);
    }
    if(((z+1)<DIMZ) && label[idx+OFFZ]==-3){//front
        label[idx+OFFZ]=-2; phi[idx+OFFZ]=-2;
        ll_pushnew(Ln2,x,y,z+1,idx+OFFZ);
    }
    if(((z-1)>=0) && label[idx-OFFZ]==-3){//back
        label[idx-OFFZ]=-2; phi[idx-OFFZ]=-2;
        ll_pushnew(Ln2,x,y,z-1,idx-OFFZ);
    }
    ll_step(Ln1);
  }

  //scan Lp1 to create Lp2
  ll_init(Lp1);
  while(Lp1->curr != NULL){
    x = Lp1->curr->x; y = Lp1->curr->y; z = Lp1->curr->z; idx = Lp1->curr->idx;
    
    if(((y+1)<DIMY) && label[idx+OFFY]==3){//up
        label[idx+OFFY]=2; phi[idx+OFFY]=2;
        ll_pushnew(Lp2,x,y+1,z,idx+OFFY);
    }
    if(((y-1)>=0) && label[idx-OFFY]==3){//down
        label[idx-OFFY]=2; phi[idx-OFFY]=2;
        ll_pushnew(Lp2,x,y-1,z,idx-OFFY);
    }
    if(((x+1)<DIMX) && label[idx+OFFX]==3){//right
        label[idx+OFFX]=2; phi[idx+OFFX]=2;
        ll_pushnew(Lp2,x+1,y,z,idx+OFFX);
    }
    if(((x-1)>=0) && label[idx-OFFX]==3){//left
        label[idx-OFFX]=2; phi[idx-OFFX]=2;
        ll_pushnew(Lp2,x-1,y,z,idx-OFFX);
    }
    if(((z+1)<DIMZ) && label[idx+OFFZ]==3){//front
        label[idx+OFFZ]=2; phi[idx+OFFZ]=2;
        ll_pushnew(Lp2,x,y,z+1,idx+OFFZ);
    }
    if(((z-1)>=0) && label[idx-OFFZ]==3){//back
        label[idx-OFFZ]=2; phi[idx-OFFZ]=2;
        ll_pushnew(Lp2,x,y,z-1,idx-OFFZ);
    }
    ll_step(Lp1);
  }
}
