#include "topology.h"

/***************************************************************************//*** 
*  \brief calculate distance to interface as well as cell marker
*******************************************************************************/
/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
 *
 * cell marker =  0: error
 *             = -1: interface in the caller cell
 *             = +1: interface in the other cell
*/

/* generic */
real Topology::distance_int(const Sign dir, const Comp & m,
                            const int i, const int j, const int k,
                            Sign & cell_marker) const {
  if        (m==Comp::i()) {
    return distance_int_x(dir,i,j,k,cell_marker);
  } else if (m==Comp::j()) {
    return distance_int_y(dir,i,j,k,cell_marker);
  } else {
    return distance_int_z(dir,i,j,k,cell_marker);
  }

  return 0.0;
}


/* x-direction */
real Topology::distance_int_x(const Sign dir, 
                              const int i, const int j, const int k,
                              Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_x(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dxc(i));
  }

  boil::aout<<"Topology::int_dist_x: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<int(dir)<<" "
            <<(*iflag)[i][j][k]<<" "<<(*iflag)[i+dir][j][k]<<" "
            <<(*clr)[i][j][k]<<" "<<(*clr)[i+dir][j][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

/* y-direction */
real Topology::distance_int_y(const Sign dir,
                              const int i, const int j, const int k,
                              Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_y(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dyc(j));
  }

  boil::aout<<"Topology::int_dist_y: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<int(dir)<<" "
            <<(*iflag)[i][j][k]<<" "<<(*iflag)[i][j+dir][k]<<" "
            <<(*clr)[i][j][k]<<" "<<(*clr)[i][j+dir][k]<<boil::endl;
  exit(0);

  return 0.0;
}

/* z-direction */
real Topology::distance_int_z(const Sign dir,
                              const int i, const int j, const int k,
                              Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_z(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dzc(k));
  }

  boil::aout<<"Topology::int_dist_z: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<int(dir)<<" "
            <<(*iflag)[i][j][k]<<" "<<(*iflag)[i][j][k+dir]<<" "
            <<(*clr)[i][j][k]<<" "<<(*clr)[i][j][k+dir]<<boil::endl;
  exit(0);

  return 0.0;
}

/******************************************************************************/
Sign Topology::distance1D_int_x(const int i, const int j, const int k,
                                const Sign dir, real & dist) const {
  real centrex = vf->xc(i);
  if(dir>0) {
    real edgex= centrex+0.5*vf->dxc(i);
    real intx = (*fs)[Comp::i()][i+1][j][k];
    real offx = vf->xc(i+1);
    if(intx>=centrex&&intx<=edgex) {
      dist = intx-centrex;
      return Sign::neg();
    } else if(intx>=edgex&&intx<=offx) {
      dist = intx-centrex;
      return Sign::pos();
    }
  } else {
    real edgex= centrex-0.5*vf->dxc(i);
    real intx = (*fs)[Comp::i()][i  ][j][k];
    real offx = vf->xc(i-1);
    if(intx<=centrex&&intx>=edgex) {
      dist = centrex-intx;
      return Sign::neg();
    } else if(intx>=offx&&intx<=edgex) {
      dist = centrex-intx;
      return Sign::pos();
    } 
  }
 
  return Sign::undefined();
}

Sign Topology::distance1D_int_y(const int i, const int j, const int k,
                                const Sign dir, real & dist) const {
  real centrey = vf->yc(j);
  if(dir>0) {
    real edgey= centrey+0.5*vf->dyc(j);
    real inty = (*fs)[Comp::j()][i][j+1][k];
    real offy = vf->yc(j+1);
    if(inty>=centrey&&inty<=edgey) {
      dist = inty-centrey;
      return Sign::neg();
    } else if(inty>=edgey&&inty<=offy) {
      dist = inty-centrey;
      return Sign::pos();
    }
  } else {
    real edgey= centrey-0.5*vf->dyc(j);
    real inty = (*fs)[Comp::j()][i][j  ][k];
    real offy = vf->yc(j-1);
    if(inty<=centrey&&inty>=edgey) {
      dist = centrey-inty;
      return Sign::neg();
    } else if(inty>=offy&&inty<=edgey) {
      dist = centrey-inty;
      return Sign::pos();
    } 
  }
 
  return Sign::undefined();
}

Sign Topology::distance1D_int_z(const int i, const int j, const int k,
                                const Sign dir, real & dist) const {
  real centrez = vf->zc(k);
  if(dir>0) {
    real edgez= centrez+0.5*vf->dzc(k);
    real intz = (*fs)[Comp::k()][i][j][k+1];
    real offz = vf->zc(k+1);
    if(intz>=centrez&&intz<=edgez) {
      dist = intz-centrez;
      return Sign::neg();
    } else if(intz>=edgez&&intz<=offz) {
      dist = intz-centrez;
      return Sign::pos();
    }
  } else {
    real edgez= centrez-0.5*vf->dzc(k);
    real intz = (*fs)[Comp::k()][i][j][k  ];
    real offz = vf->zc(k-1);
    if(intz<=centrez&&intz>=edgez) {
      dist = centrez-intz;
      return Sign::neg();
    } else if(intz>=offz&&intz<=edgez) {
      dist = centrez-intz;
      return Sign::pos();
    } 
  }
 
  return Sign::undefined();
}
