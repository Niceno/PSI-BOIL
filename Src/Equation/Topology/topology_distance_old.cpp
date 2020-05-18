#include "topology.h"

/***************************************************************************//*** 
*  \brief calculate distance to interface as well as cell marker using old vals
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
real Topology::distance_int_old(const Sign dir, const Comp & m,
                                const int i, const int j, const int k,
                                Sign & cell_marker) const {
  if        (m==Comp::i()) {
    return distance_int_x_old(dir,i,j,k,cell_marker);
  } else if (m==Comp::j()) {
    return distance_int_y_old(dir,i,j,k,cell_marker);
  } else {
    return distance_int_z_old(dir,i,j,k,cell_marker);
  }

  return 0.0;
}


/* x-direction */
real Topology::distance_int_x_old(const Sign dir, 
                                  const int i, const int j, const int k,
                                  Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_x_old(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dxc(i));
  }

  boil::aout<<"Topology::int_dist_x_old: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflagold[i][j][k]<<" "<<iflagold[i+dir][j][k]<<" "
            <<clrold[i][j][k]<<" "<<clrold[i+dir][j][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

/* y-direction */
real Topology::distance_int_y_old(const Sign dir,
                                  const int i, const int j, const int k,
                                  Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_y_old(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dyc(j));
  }

  boil::aout<<"Topology::int_dist_y_old: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflagold[i][j][k]<<" "<<iflagold[i][j+dir][k]<<" "
            <<clrold[i][j][k]<<" "<<clrold[i][j+dir][k]<<boil::endl;
  exit(0);

  return 0.0;
}

/* z-direction */
real Topology::distance_int_z_old(const Sign dir,
                                  const int i, const int j, const int k,
                                  Sign & cell_marker) const {
  real dist;

  cell_marker = distance1D_int_z_old(i,j,k,dir,dist);
  if(cell_marker != Sign::undefined()) {
    return std::max(dist,close_to_cc*clrold.dzc(k));
  }

  boil::aout<<"Topology::int_dist_z_old: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflagold[i][j][k]<<" "<<iflagold[i][j][k+dir]<<" "
            <<clrold[i][j][k]<<" "<<clrold[i][j][k+dir]<<boil::endl;
  exit(0);

  return 0.0;
}

/******************************************************************************/
Sign Topology::distance1D_int_x_old(const int i, const int j, const int k,
                                    const Sign dir, real & dist) const {
  real centrex = vf->xc(i);
  if(dir>0) {
    real edgex= centrex+0.5*vf->dxc(i);
    real intx = fsold[Comp::i()][i+1][j][k];
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
    real intx = fsold[Comp::i()][i  ][j][k];
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

Sign Topology::distance1D_int_y_old(const int i, const int j, const int k,
                                    const Sign dir, real & dist) const {
  real centrey = vf->yc(j);
  if(dir>0) {
    real edgey= centrey+0.5*vf->dyc(j);
    real inty = fsold[Comp::j()][i][j+1][k];
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
    real inty = fsold[Comp::j()][i][j  ][k];
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

Sign Topology::distance1D_int_z_old(const int i, const int j, const int k,
                                    const Sign dir, real & dist) const {
  real centrez = vf->zc(k);
  if(dir>0) {
    real edgez= centrez+0.5*vf->dzc(k);
    real intz = fsold[Comp::k()][i][j][k+1];
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
    real intz = fsold[Comp::k()][i][j][k  ];
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
