#ifndef UTILS_H
#define UTILS_H

#ifndef DEBUG
#define DEBUG 2
#endif
#include "mex.h"
#include "linkedList.cpp"
#include "array.h"
#include "helperMEX.h"
#include "debugMEX.h"
// --1. mex I/O functions
inline
mxArray *mxReadField(const mxArray *mstruct,const char *fieldstr) {
   mxArray *result;
   result = mxGetField(mstruct,0,fieldstr);
   if (result == NULL) {
      mexPrintf(fieldstr);
      mexErrMsgTxt(" field missing.");
   }
   if (DEBUG>=3) mexPrintf("Read %s.\n",fieldstr);
   return result;
}

inline
void mxWriteField(mxArray *mstruct,const char *fieldstr,mxArray *mfield) {
   if ( mxGetField(mstruct,0,fieldstr) != NULL )
      mxDestroyArray(mxGetField(mstruct,0,fieldstr));
   mxSetField(mstruct,0,fieldstr,mfield);
   if (DEBUG>=3) mexPrintf("Write %s.\n",fieldstr);
}

inline
void mxWriteField(mxArray *mstruct,int ind,const char *fieldstr,mxArray *mfield) {
   if ( mxGetField(mstruct,ind,fieldstr) != NULL )
      mxDestroyArray(mxGetField(mstruct,ind,fieldstr));
   mxSetField(mstruct,ind,fieldstr,mfield);
   if (DEBUG>=3) mexPrintf("Write %s.\n",fieldstr);
}

inline
double mxReadScalar(const mxArray *mscalar) {
   if (DEBUG>=3) mexPrintf("Value = %g.\n",*mxGetPr(mscalar));
   if (mscalar==NULL || !mxIsDouble(mscalar)){
      mexErrMsgTxt("Sth. wrong with the Scalar\n");
      return 0;
   }else{

      return (*mxGetPr(mscalar));
   }
}

inline
mxArray *mxWriteScalar(double var) {
   mxArray *result;
   if (DEBUG>=3) mexPrintf("Value = %g.\n",var);
   result = mxCreateDoubleMatrix(1,1,mxREAL);
   *mxGetPr(result) = var;
   return result;
}

inline
mxArray *mxWriteScalar(unsigned long var) {
   mxArray *result;
   if (DEBUG>=3) mexPrintf("Value = %g.\n",var);
   result = mxCreateNumericMatrix(1,1,mxUINT64_CLASS,mxREAL);
   *((unsigned long*)(mxGetPr(result))) = var;
   return result;
}

inline
mxArray *mxWriteScalar(bool var) {
   mxArray *result;
   if (DEBUG>=3) mexPrintf("Value = %g.\n",var);
   result = mxCreateLogicalMatrix(1,1);
   *((bool*)(mxGetPr(result))) = var;
   return result;
}

template <typename T>
T* mxReadVector(const mxArray *mvector) {
   double *mdouble;
   T *result;
   int num = mxGetNumberOfElements(mvector);
   result = new T[num];
   mdouble = mxGetPr(mvector);
   for (int ii = 0 ; ii < num ; ii++ ){
      result[ii] = mdouble[ii];
   }
   if (DEBUG>=3) {
      mexPrintf("Value = ");
      for (int ii = 0 ; ii < num ; ii++ )
      //mexPrintf(str,result[ii]); //? how to print it
      mexPrintf("\n");
   }
   return result;
}
template <typename T>
void mxReadVector2(const mxArray *mvector,T *result) {
   double *mdouble;
   int num = mxGetNumberOfElements(mvector);
   memcpy(result,mxGetPr(mvector),sizeof(T)*num);
   if (DEBUG>=3) {
      mexPrintf("Value = ");
      for (int ii = 0 ; ii < num ; ii++ )
      mexPrintf("%g ",result[ii]);
      mexPrintf("\n");
   }
}
// dynamic array
template <typename T>
mxArray *mxWriteVector(int mm,int nn,T *var) {
   mxArray *result;
   double *mdouble;
   result = mxCreateDoubleMatrix(mm,nn,mxREAL);
   mdouble = mxGetPr(result);
   for (int ii = 0 ; ii < mm*nn ; ii++ ){mdouble[ii] = var[ii];}
   mxFree(var);
   return result;
}


// constant array
template <typename T>
mxArray *mxWriteVector2(int mm,int nn,T *var) {
   mxArray *result;
   double *mdouble;
   result = mxCreateDoubleMatrix(mm,nn,mxREAL);
   mdouble = mxGetPr(result);
   for (int ii = 0 ; ii < mm*nn ; ii++ ){mdouble[ii] = var[ii];}
   return result;
}

// dynamic array
template <typename T>
mxArray *mxWriteVectorPtr(int mm,int nn,T *var) {
   mxArray *result;
   double *mdouble;
   result = mxCreateDoubleMatrix(mm,nn,mxREAL);
   mdouble = mxGetPr(result);
   for (int ii = 0 ; ii < mm*nn ; ii++ ){mdouble[ii] = var[ii];}
   delete[] var;
   return result;
}


//  --2. Handy Algorithms
inline
int U_randmult(double *weight, int len) {
	double sum = 0.0, mass;
	int cc = 0;
	for (int i =0 ; i < len ; i ++)
		sum += weight[i];
	mass = drand48() * sum;
	//printf("ss:%f,%f\n",mass,sum);
	while (1) {
		mass -= weight[cc];
		if ( mass <= 0.0 ) break;
		cc ++;
	}
	return cc;
}

template <typename T>
void U_swap(T* arr,int x,int y) {
	T *t1,*t2,tt;
	t1 = arr+x;
	t2 = arr+y;
	tt = *t1;
	*t1 = *t2;
	*t2 = tt;
}
	//ascending order
	// -- return sorted array: *list
	// -- corresponding index of original position: *ind
template <typename T>
void U_quicksort(T *list,int *ind,int m,int n) {
	int i,j,k;
	double key;
	if( m < n) {
		k = (int)(m+n)/2;
		U_swap<T>(list,m,k);
		U_swap<int>(ind,m,k);
		key = list[m];
		i = m+1;
		j = n;
		while(i <= j) {
			while((i <= n) && (list[i] <= key))
				i++;
			while((j >= m) && (list[j] > key))
				j--;
			if( i < j) {
				U_swap<T>(list,i,j);
				U_swap<int>(ind,i,j);
			}
		}
		// swap two elements
		U_swap<T>(list,m,j);
		U_swap<int>(ind,m,j);
		// recursively sort the lesser list
		U_quicksort<T>(list,ind,m,j-1);
		U_quicksort<T>(list,ind,j+1,n);
	}
}
inline
int U_FindArray(int index, linkedList<int> &check_labels){
	int found = -1,pos=0;
    linkedListNode<int> *node;
    node = check_labels.getFirst();
    if(node==NULL)
        mexErrMsgTxt("FIndArray: check_labels is empty...\n");
    while(node != NULL){
		if(node->getData() == index){
			found = pos;
			break;
		}
        pos++;
        node = node->getNext();
	}
	return found;
}
inline
int U_FindArray(int index, int* check_labels, int num_check){
	int found = -1;
	for( int j = 0; j < num_check; j++ ){
		if(check_labels[j] == index){
			found = j;
			break;
		}
	}
	return found;
}
template <typename T>
void PrintList(linkedList<T> ll) {
    linkedListNode<T> * node= ll.getFirst();
         while(node!=NULL){
             std::cout<<node->getData()<<",";
         node = node->getNext();
         }
         std::cout<<std::endl;
}
#endif
