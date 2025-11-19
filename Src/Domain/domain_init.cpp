#include "domain.h"
#include "../Parallel/communicator.h"
#include "../Plot/plot.h"

/******************************************************************************/
void Domain::init(const Decompose & dec) {

  const int ndims = 3;

  dims       = new int[3];
  coords     = new int[3];
  neighbours = new int[6];
	
  /* initialize dims, coords and neighbours */
  for(int i=0; i<ndims; i++) dims      [i] = 1;
  for(int i=0; i<ndims; i++) coords    [i] = 0;
  for(int i=0; i<6;     i++) neighbours[i] = par_proc_null;

  if( dec == Decompose::no() ) return;

  /* get the resolution in each direction */
  //int res[] = { gi()-2 , gj()-2 , gk()-2 };
  int res[] = { gi()-2*boil::BW , gj()-2*boil::BW , gk()-2*boil::BW }; //2019.0615

  /* take care of constrained decompositions */
  if( dec != Decompose::x()  && 
      dec != Decompose::xy() && 
      dec != Decompose::yx() && 
      dec != Decompose::xz() && 
      dec != Decompose::zx() && 
      dec != Decompose::xyz() ) res[0] = 1;
  if( dec != Decompose::y()  && 
      dec != Decompose::yx() && 
      dec != Decompose::xy() && 
      dec != Decompose::yz() && 
      dec != Decompose::zy() && 
      dec != Decompose::xyz() ) res[1] = 1;
  if( dec != Decompose::z()  && 
      dec != Decompose::xz() && 
      dec != Decompose::zx() && 
      dec != Decompose::yz() && 
      dec != Decompose::zy() && 
      dec != Decompose::xyz() ) res[2] = 1;

  if(boil::cart.nproc()>1)
    distribute(boil::cart.nproc(), ndims, dims, res);

  /* debugging info 
  boil::aout << boil::pid << " dims[i]: ";
  for(int i=0; i<ndims; i++)
    boil::aout << " " << dims[i];
  boil::aout << boil::endl;
  */

  /* find the neighbours */
  for(int i=0; i<dims[0]; i++)
    for(int j=0; j<dims[1]; j++)
      for(int k=0; k<dims[2]; k++) {
        int tag = i*dims[1]*dims[2] + j*dims[2] + k;
    
        if(tag == boil::cart.iam()) {

          if(i>0        ) neighbours[Dir::imin()] = tag - dims[1]*dims[2];
          if(i<dims[0]-1) neighbours[Dir::imax()] = tag + dims[1]*dims[2];
          if(j>0        ) neighbours[Dir::jmin()] = tag - dims[2];
          if(j<dims[1]-1) neighbours[Dir::jmax()] = tag + dims[2];
          if(k>0        ) neighbours[Dir::kmin()] = tag - 1;      
          if(k<dims[2]-1) neighbours[Dir::kmax()] = tag + 1; 

          if(period(0)) {
            if(i==0) 
              neighbours[Dir::imin()] = (dims[0]-1)*dims[1]*dims[2] + j*dims[2] + k;
            if(i==dims[0]-1) 
              neighbours[Dir::imax()] = j*dims[2] + k;
          }
          if(period(1)) {
            if(j==0) 
              neighbours[Dir::jmin()] = i*dims[1]*dims[2] + (dims[1]-1)*dims[2] + k;
            if(j==dims[1]-1) 
              neighbours[Dir::jmax()] = i*dims[1]*dims[2] + k;
          }
          if(period(2)) {
            if(k==0) 
              neighbours[Dir::kmin()] = i*dims[1]*dims[2] + j*dims[2] + dims[2]-1;
            if(k==dims[2]-1) 
              neighbours[Dir::kmax()] = i*dims[1]*dims[2] + j*dims[2];
          }
        }
      }

  /* debugging info 
  boil::aout << boil::pid << " neighbours: ";
  for(int i=0; i<6; i++) 
    boil::aout << " " << neighbours[i];
  boil::aout << boil::endl;
  */
  
  /* find your coordinates */
  for(int i=0; i<dims[0]; i++)
    for(int j=0; j<dims[1]; j++)
      for(int k=0; k<dims[2]; k++) {
        int tag = i*dims[1]*dims[2] + j*dims[2] + k;
        if(tag == boil::cart.iam()) {
          coords[0] = i;
          coords[1] = j;
          coords[2] = k;
        }
      }

  /* debugging info 
  boil::aout << boil::pid << " coords[i]: ";
  for(int i=0; i<ndims; i++)
    boil::aout << " " << coords[i];
  boil::aout << boil::endl;
  */
}	
