#include "phasechange4.h"

/******************************************************************************/
real PhaseChange4::distance_center(const Sign sig, const Comp & m,
                                 const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell center  
*******************************************************************************/
  if       (m==Comp::i()) {
    if(sig<0) {
      return phi.dxw(i);
    } else {
      return phi.dxe(i);
    }
  } else if(m==Comp::j()) {
    if(sig<0) {
      return phi.dys(j);
    } else {
      return phi.dyn(j);
    }
  } else {
    if(sig<0) {
      return phi.dzb(k);
    } else {
      return phi.dzt(k);
    }
  }

  return 0.0;
}

/******************************************************************************/
real PhaseChange4::distance_face(const Sign sig, const Comp & m,
                                 const int i, const int j, const int k) {
/***************************************************************************//*** 
*  \brief calculate distance to neighboring cell face  
*******************************************************************************/
  if       (m==Comp::i()) {
    return 0.5*phi.dxc(i);
  } else if(m==Comp::j()) {
    return 0.5*phi.dyc(j);
  } else {
    return 0.5*phi.dzc(k);
  }

  return 0.0;
}

/***************************************************************************//*** 
*  \brief calculate distance to interface, as well as interface tpr  
*******************************************************************************/
/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
*/

/* generic */
real PhaseChange4::distance_int(const Sign dir, const Comp & m,
                                const int i, const int j, const int k,
                                real & tint) {
  if        (m==Comp::i()) {
    return distance_int_x(dir,i,j,k,tint);
  } else if (m==Comp::j()) {
    return distance_int_y(dir,i,j,k,tint);
  } else {
    return distance_int_z(dir,i,j,k,tint);
  }

  return 0.0;
}


/* x-direction */
real PhaseChange4::distance_int_x(const Sign dir, 
                                  const int i, const int j, const int k,
                                  real & tint) {
  real dist;

  if(distance1D_int_x(i,j,k,dir,tint,dist)) {
    return dist;
  }

  boil::aout<<"PhaseChange4::int_dist_x: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i+dir][j][k]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i+dir][j][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChange4::distance1D_int_x(const int i, const int j, const int k,
                                    const Sign dir, real & tint, real & dist) {
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
real PhaseChange4::distance_int_y(const Sign dir, 
                                  const int i, const int j, const int k,
                                  real & tint) {
  real dist;

  if(distance1D_int_y(i,j,k,dir,tint,dist)) {
    return dist;
  }

  boil::aout<<"PhaseChange4::int_dist_y: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i][j+dir][k]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i][j+dir][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChange4::distance1D_int_y(const int i, const int j, const int k,
                                    const Sign dir, real & tint, real & dist) {
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
real PhaseChange4::distance_int_z(const Sign dir, 
                                  const int i, const int j, const int k,
                                  real & tint) {
  real dist;

  if(distance1D_int_z(i,j,k,dir,tint,dist)) {
    return dist;
  }

  boil::aout<<"PhaseChange4::int_dist_z: Error! Flag inconsistent w/ vol. fraction!\n";
  boil::aout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
            <<iflag[i][j][k]<<" "<<iflag[i][j][k+dir]<<" "
            <<(clr)[i][j][k]<<" "<<(clr)[i][j][k+dir]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool PhaseChange4::distance1D_int_z(const int i, const int j, const int k,
                                    const Sign dir, real & tint, real & dist) {
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
