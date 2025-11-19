#include <cstdlib>

#include "global_precision.h"

void alloc1d(real   ** v, const int n); 
void alloc2d(real  *** v, const int ni, const int nj); 
void alloc3d(real **** v, const int ni, const int nj, const int nk);

void alloc1d(int   ** v, const int n); 
void alloc2d(int  *** v, const int ni, const int nj); 
void alloc3d(int **** v, const int ni, const int nj, const int nk);

void alloc1d(bool   ** v, const int n); 
void alloc2d(bool  *** v, const int ni, const int nj); 
void alloc3d(bool **** v, const int ni, const int nj, const int nk);

void dealloc1d(real   ** v); 
void dealloc2d(real  *** v); 
void dealloc3d(real **** v);

void dealloc1d(int   ** v); 
void dealloc2d(int  *** v); 
void dealloc3d(int **** v);

void dealloc1d(bool   ** v); 
void dealloc2d(bool  *** v); 
void dealloc3d(bool **** v);
