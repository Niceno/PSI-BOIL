#include "phasechangevof.h"

/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
*/

/* x-direction */
real PhaseChangeVOF::distance_x(const int i, const int j, const int k,
                               const int dir, real & tint) {
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(distance1D_x(i,j,k,of,tint,dist)) {
    return dist;
  }

  boil::oout<<"PhaseChangeVOF::int_dist_x: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i+of][j][k]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i+of][j][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChangeVOF::distance1D_x(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrex = phi.xc(i);
  if(dir>0) {
    real edgex= centrex+0.5*phi.dxc(i);
    real intx = (fs)[Comp::i()][i+1][j][k];
    real offx = phi.xc(i+1);
    if(intx>=centrex&&intx<=edgex) {
      tint = Tint(i,j,k);
      dist = intx-centrex;
      return true;
    } else if(intx>=edgex&&intx<=offx) {
      tint = Tint(i+1,j,k);
      dist = intx-centrex;
      return true;
    }
  } else {
    real edgex= centrex-0.5*phi.dxc(i);
    real intx = (fs)[Comp::i()][i  ][j][k];
    real offx = phi.xc(i-1);
    if(intx<=centrex&&intx>=edgex) {
      tint = Tint(i,j,k);
      dist = centrex-intx;
      return true;
    } else if(intx>=offx&&intx<=edgex) {
      tint = Tint(i-1,j,k);
      dist = centrex-intx;
      return true;
    } 
  }
 
  return false;
}

/* y-direction */
real PhaseChangeVOF::distance_y(const int i, const int j, const int k,
                               const int dir, real & tint) {
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(distance1D_y(i,j,k,of,tint,dist)) {
    return dist;
  }

  boil::oout<<"PhaseChangeVOF::int_dist_y: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i][j+of][k]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i][j+of][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChangeVOF::distance1D_y(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrey = phi.yc(j);
  if(dir>0) {
    real edgey= centrey+0.5*phi.dyc(j);
    real inty = (fs)[Comp::j()][i][j+1][k];
    real offy = phi.yc(j+1);
    if(inty>=centrey&&inty<=edgey) {
      tint = Tint(i,j,k);
      dist = inty-centrey;
      return true;
    } else if(inty>=edgey&&inty<=offy) {
      tint = Tint(i,j+1,k);
      dist = inty-centrey;
      return true;
    }
  } else {
    real edgey= centrey-0.5*phi.dyc(j);
    real inty = (fs)[Comp::j()][i][j  ][k];
    real offy = phi.yc(j-1);
    if(inty<=centrey&&inty>=edgey) {
      tint = Tint(i,j,k);
      dist = centrey-inty;
      return true;
    } else if(inty>=offy&&inty<=edgey) {
      tint = Tint(i,j-1,k);
      dist = centrey-inty;
      return true;
    } 
  }
 
  return false;
}

/* z-direction */
real PhaseChangeVOF::distance_z(const int i, const int j, const int k,
                                const int dir, real & tint) {
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(distance1D_z(i,j,k,of,tint,dist)) {
    return dist;
  }

  boil::oout<<"PhaseChangeVOF::int_dist_z: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i][j][k+of]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i][j][k+of]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChangeVOF::distance1D_z(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrez = phi.zc(k);
  if(dir>0) {
    real edgez= centrez+0.5*phi.dzc(k);
    real intz = (fs)[Comp::k()][i][j][k+1];
    real offz = phi.zc(k+1);
    if(intz>=centrez&&intz<=edgez) {
      tint = Tint(i,j,k);
      dist = intz-centrez;
      return true;
    } else if(intz>=edgez&&intz<=offz) {
      tint = Tint(i,j,k+1);
      dist = intz-centrez;
      return true;
    }
  } else {
    real edgez= centrez-0.5*phi.dzc(k);
    real intz = (fs)[Comp::k()][i][j][k  ];
    real offz = phi.zc(k-1);
    if(intz<=centrez&&intz>=edgez) {
      tint = Tint(i,j,k);
      dist = centrez-intz;
      return true;
    } else if(intz>=offz&&intz<=edgez) {
      tint = Tint(i,j,k-1);
      dist = centrez-intz;
      return true;
    } 
  }
 
  return false;
}
