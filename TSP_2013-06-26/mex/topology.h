#ifndef TOPOLOGY_H
#define TOPOLOGY_H

inline int calc_T_table(bool* neighborhood, double* table)
{
   int index = 0;
   for (int i=7; i>=0; i--)
   {
      index <<= 1;
      if (neighborhood[i])
         index++;
   }
   /*if (table[index]<0)
      mexErrMsgTxt("Topological Table Error");

   mexPrintf("\nNEIGHBORHOOD: %d  ", index);
   for (int i=0; i<8; i++)
      mexPrintf("%d ", neighborhood[i]);
   mexPrintf("\n%d\n", (int)(table[index]));*/
   return (int)(table[index]);
}
inline int topology_number_label(int index, int check_label, int* label, bool* neighborhood, int xdim, int ydim, double* T4Table)
{
   int x = index%xdim;
   int y = index/xdim;

   if (check_label<0) // dummy
      return 1;

   neighborhood[1] = (y>0 && label[index-xdim]==check_label);
   neighborhood[3] = (x<xdim-1 && label[index+1]==check_label);
   neighborhood[5] = (y<ydim-1 && label[index+xdim]==check_label);
   neighborhood[7] = (x>0 && label[index-1]==check_label);
   neighborhood[0] = ((neighborhood[1] || neighborhood[7]) && (x>0 && y>0 && label[index-1-xdim]==check_label));
   neighborhood[2] = ((neighborhood[1] || neighborhood[3]) && (x<xdim-1 && y>0 && label[index+1-xdim]==check_label));
   neighborhood[4] = ((neighborhood[3] || neighborhood[5]) && (x<xdim-1 && y<ydim-1 && label[index+1+xdim]==check_label));
   neighborhood[6] = ((neighborhood[5] || neighborhood[7]) && (x>0 && y<ydim-1 && label[index-1+xdim]==check_label));
   return (int)(calc_T_table(neighborhood, T4Table));
}
inline bool check_topology(int index, int* label, bool* neighborhood, int xdim, int ydim, double* T4Table)
{
   int x = index%xdim;
   int y = index/xdim;

   // check the left label
   bool top_ok =
      ( x<=0 || (topology_number_label(index, label[index-1], label, neighborhood, xdim, ydim, T4Table)==1) ) &&
      ( y<=0 || (topology_number_label(index, label[index-xdim], label, neighborhood, xdim, ydim, T4Table)==1) ) &&
      ( x>=xdim-1 || (topology_number_label(index, label[index+1], label, neighborhood, xdim, ydim, T4Table)==1) ) &&
      ( y>=ydim-1 || (topology_number_label(index, label[index+xdim], label, neighborhood, xdim, ydim, T4Table)==1) );

   return top_ok;


   int tmp=topology_number_label(index, label[index], label, neighborhood, xdim, ydim, T4Table);
   tmp=0;
   for(int i=0;i<8;i++){tmp+= neighborhood[i];}
   if(tmp==1){
      //if(neighborhood[7]==1&&label[index]==87&&y==70){printf("gotttaaa,%d,%d,%d\n",x,y,index);}
        //return true;
        return top_ok;
      }else{
   return top_ok;}
}

#endif
