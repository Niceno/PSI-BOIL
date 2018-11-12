#include "global_malloc.h"

/******************************************************************************/
void alloc1d(real ** v, const int n) {
/*----------------------------------------+
|  allocate memory as a contiguous block  |
+----------------------------------------*/	
	
  (*v) = new real [n];

  /* initialize to zero: very important for MPI version */
  for(int i=0; i<n; i++)
    (*v)[i] = 0;
}	

/******************************************************************************/
void alloc2d(real *** v, const int ni, const int nj) {
/*--------------------------------------------------------+
|  allocate memory for a 2D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new real * [ni];
  
  (*v)[0] = new real [ni * nj]; // important: contiguous block
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = 0.0;
}	

/******************************************************************************/
void alloc3d(real **** v, const int ni, const int nj, const int nk) {
/*--------------------------------------------------------+
|  allocate memory for a 3D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new real ** [ni];
  
  (*v)[0] = new real *  [ni * nj];
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;

  (*v)[0][0] = new real [ni * nj * nk]; // important: contiguous block
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = (*v)[0][0] + i*nj*nk + j*nk;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      for(int k=0; k<nk; k++)
        (*v)[i][j][k] = 0.0;
}

/******************************************************************************/
void alloc1d(int ** v, const int n) {
/*----------------------------------------+
|  allocate memory as a contiguous block  |
+----------------------------------------*/	
	
  (*v) = new int [n];

  /* initialize to zero: very important for MPI version */
  for(int i=0; i<n; i++)
    (*v)[i] = 0;
}	

/******************************************************************************/
void alloc2d(int *** v, const int ni, const int nj) {
/*--------------------------------------------------------+
|  allocate memory for a 2D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new int * [ni];
  
  (*v)[0] = new int [ni * nj]; // important: contiguous block
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = 0;
}	

/******************************************************************************/
void alloc3d(int **** v, const int ni, const int nj, const int nk) {
/*--------------------------------------------------------+
|  allocate memory for a 3D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new int ** [ni];
  
  (*v)[0] = new int *  [ni * nj];
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;

  (*v)[0][0] = new int [ni * nj * nk]; // important: contiguous block
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = (*v)[0][0] + i*nj*nk + j*nk;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      for(int k=0; k<nk; k++)
        (*v)[i][j][k] = 0;
}

/******************************************************************************/
void alloc1d(bool ** v, const int n) {
/*----------------------------------------+
|  allocate memory as a contiguous block  |
+----------------------------------------*/	
	
  (*v) = new bool [n];

  /* initialize to zero: very important for MPI version */
  for(int i=0; i<n; i++)
    (*v)[i] = 0;
}	

/******************************************************************************/
void alloc2d(bool *** v, const int ni, const int nj) {
/*--------------------------------------------------------+
|  allocate memory for a 2D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new bool * [ni];
  
  (*v)[0] = new bool [ni * nj]; // important: contiguous block
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = 0;
}	

/******************************************************************************/
void alloc3d(bool **** v, const int ni, const int nj, const int nk) {
/*--------------------------------------------------------+
|  allocate memory for a 3D array, as a contiguous block  |
+--------------------------------------------------------*/	
	
  (*v) = new bool ** [ni];
  
  (*v)[0] = new bool *  [ni * nj];
  for(int i=0; i<ni; i++)
    (*v)[i] = (*v)[0] + i*nj;

  (*v)[0][0] = new bool [ni * nj * nk]; // important: contiguous block
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      (*v)[i][j] = (*v)[0][0] + i*nj*nk + j*nk;
  
  /* initialize to zero: very important for MPI version */
  for(int i=0; i<ni; i++)
    for(int j=0; j<nj; j++)
      for(int k=0; k<nk; k++)
        (*v)[i][j][k] = 0;
}

/******************************************************************************/
void dealloc1d(real ** v) {
  delete [] (*v); (*v) = NULL;
}

/******************************************************************************/
void dealloc2d(real *** v) {
  delete [] (*v)[0]; (*v)[0] = NULL;
  delete [] (*v);    (*v)    = NULL;
}

/******************************************************************************/
void dealloc3d(real **** v) {
  delete [] (*v)[0][0]; (*v)[0][0] = NULL;
  delete [] (*v)[0];    (*v)[0]    = NULL;
  delete [] (*v);       (*v)       = NULL;
}

/******************************************************************************/
void dealloc1d(int ** v) {
  delete [] (*v); (*v) = NULL;
}

/******************************************************************************/
void dealloc2d(int *** v) {
  delete [] (*v)[0]; (*v)[0] = NULL;
  delete [] (*v);    (*v)    = NULL;
}

/******************************************************************************/
void dealloc3d(int **** v) {
  delete [] (*v)[0][0]; (*v)[0][0] = NULL;
  delete [] (*v)[0];    (*v)[0]    = NULL;
  delete [] (*v);       (*v)       = NULL;
}

/******************************************************************************/
void dealloc1d(bool ** v) {
  delete [] (*v); (*v) = NULL;
}

/******************************************************************************/
void dealloc2d(bool *** v) {
  delete [] (*v)[0]; (*v)[0] = NULL;
  delete [] (*v);    (*v)    = NULL;
}

/******************************************************************************/
void dealloc3d(bool **** v) {
  delete [] (*v)[0][0]; (*v)[0][0] = NULL;
  delete [] (*v)[0];    (*v)[0]    = NULL;
  delete [] (*v);       (*v)       = NULL;
}
