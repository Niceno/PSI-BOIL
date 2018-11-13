#include "enthalpytif.h"

/* 
 * dir > 0: positive direction
 * dir < 0: negative direction
*/

/* x-direction */
real EnthalpyTIF::distance_x(const int i, const int j, const int k,
                             const int dir, real & tint,
                             const bool old) { 
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(!old) {
    if(distance1D_x(i,j,k,of,tint,dist)) {
        //if(j==1&&k==1) boil::oout<<"ETIF::1Dx "<<i<<" "<<of<<" "<<dist<<" "<<(*clr)[i][j][k]<<" "<<(*clr)[i+of][j][k]<<" "<<(*fs)[Comp::i()][i+(of>0)][j][k]<<boil::endl;
        return dist;
    }
  /* old? */
  } else {
    if(distance1D_xold(i,j,k,of,tint,dist)) {
      return dist;
    }
  }

  boil::oout<<"EnthalpyTIF::int_dist_x: Error! Flag inconsistent w/ vol. fraction!\n";
  if(!old)
    boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              //<<intflag[i][j][k]<<" "<<intflag[i+of][j][k]<<" "
              <<(*clr)[i][j][k]<<" "<<(*clr)[i+of][j][k]<<boil::endl;
  else
    boil::oout<<"old "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              //<<intflagold[i][j][k]<<" "<<intflagold[i+of][j][k]<<" "
              <<clrold[i][j][k]<<" "<<clrold[i+of][j][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool EnthalpyTIF::distance1D_x(const int i, const int j, const int k,
                               const int dir, real & tint, real & dist) {
  real centrex = phi.xc(i);
  if(dir>0) {
    real edgex= centrex+0.5*phi.dxc(i);
    real intx = (*fs)[Comp::i()][i+1][j][k];
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
    real intx = (*fs)[Comp::i()][i  ][j][k];
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

bool EnthalpyTIF::distance1D_xold(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrex = phi.xc(i);
  if(dir>0) {
    real edgex= centrex+0.5*phi.dxc(i);
    real intx = fsold[Comp::i()][i+1][j][k];
    real offx = phi.xc(i+1);
    if(intx>=centrex&&intx<=edgex) {
      tint = Tint_old(i,j,k);
      dist = intx-centrex;
      return true;
    } else if(intx>=edgex&&intx<=offx) {
      tint = Tint_old(i+1,j,k);
      dist = intx-centrex;
      return true;
    }
  } else {
    real edgex= centrex-0.5*phi.dxc(i);
    real intx = fsold[Comp::i()][i  ][j][k];
    real offx = phi.xc(i-1);
    if(intx<=centrex&&intx>=edgex) {
      tint = Tint_old(i,j,k);
      dist = centrex-intx;
      return true;
    } else if(intx>=offx&&intx<=edgex) {
      tint = Tint_old(i-1,j,k);
      dist = centrex-intx;
      return true;
    } 
  }
 
  return false;
}

/* y-direction */
real EnthalpyTIF::distance_y(const int i, const int j, const int k,
                             const int dir, real & tint,
                             const bool old) { 
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(!old) {
    if(distance1D_y(i,j,k,of,tint,dist)) {
        return dist;
    }
  /* old? */
  } else {
    if(distance1D_yold(i,j,k,of,tint,dist)) {
      return dist;
    }
  }

  boil::oout<<"EnthalpyTIF::int_dist_y: Error! Flag inconsistent w/ vol. fraction!\n";
  if(!old)
    boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              <<phi.yc(j)<<" "<<phi.yc(j+of)<<" "<<(*fs)[Comp::j()][i][j+(of>0)][k]<<" "
              //<<intflag[i][j][k]<<" "<<intflag[i][j+of][k]<<" "
              <<(*clr)[i][j][k]<<" "<<(*clr)[i][j+of][k]<<boil::endl;
  else
    boil::oout<<"old "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              //<<intflagold[i][j][k]<<" "<<intflagold[i][j+of][k]<<" "
              <<clrold[i][j][k]<<" "<<clrold[i][j+of][k]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool EnthalpyTIF::distance1D_y(const int i, const int j, const int k,
                               const int dir, real & tint, real & dist) {
  real centrey = phi.yc(j);
  if(dir>0) {
    real edgey= centrey+0.5*phi.dyc(j);
    real inty = (*fs)[Comp::j()][i][j+1][k];
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
    real inty = (*fs)[Comp::j()][i][j  ][k];
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

bool EnthalpyTIF::distance1D_yold(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrey = phi.yc(j);
  if(dir>0) {
    real edgey= centrey+0.5*phi.dyc(j);
    real inty = fsold[Comp::j()][i][j+1][k];
    real offy = phi.yc(j+1);
    if(inty>=centrey&&inty<=edgey) {
      tint = Tint_old(i,j,k);
      dist = inty-centrey;
      return true;
    } else if(inty>=edgey&&inty<=offy) {
      tint = Tint_old(i,j+1,k);
      dist = inty-centrey;
      return true;
    }
  } else {
    real edgey= centrey-0.5*phi.dyc(j);
    real inty = fsold[Comp::j()][i][j  ][k];
    real offy = phi.yc(j-1);
    if(inty<=centrey&&inty>=edgey) {
      tint = Tint_old(i,j,k);
      dist = centrey-inty;
      return true;
    } else if(inty>=offy&&inty<=edgey) {
      tint = Tint_old(i,j-1,k);
      dist = centrey-inty;
      return true;
    } 
  }
 
  return false;
}

/* z-direction */
real EnthalpyTIF::distance_z(const int i, const int j, const int k,
                             const int dir, real & tint,
                             const bool old) { 
  int of(1);
  if(dir<0) of = -1;

  real dist;

  if(!old) {
    if(distance1D_z(i,j,k,of,tint,dist)) {
        return dist;
    }
  /* old? */
  } else {
    if(distance1D_zold(i,j,k,of,tint,dist)) {
      return dist;
    }
  }

  boil::oout<<"EnthalpyTIF::int_dist_z: Error! Flag inconsistent w/ vol. fraction!\n";
  if(!old)
    boil::oout<<"new "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              //<<intflag[i][j][k]<<" "<<intflag[i][j][k+of]<<" "
              <<(*clr)[i][j][k]<<" "<<(*clr)[i][j][k+of]<<boil::endl;
  else
    boil::oout<<"old "<<i<<" "<<j<<" "<<k<<" "<<dir<<" "
              //<<intflagold[i][j][k]<<" "<<intflagold[i][j][k+of]<<" "
              <<clrold[i][j][k]<<" "<<clrold[i][j][k+of]<<boil::endl;
  exit(0);
 
  return 0.0;
}

bool EnthalpyTIF::distance1D_z(const int i, const int j, const int k,
                               const int dir, real & tint, real & dist) {
  real centrez = phi.zc(k);
  if(dir>0) {
    real edgez= centrez+0.5*phi.dzc(k);
    real intz = (*fs)[Comp::k()][i][j][k+1];
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
    real intz = (*fs)[Comp::k()][i][j][k  ];
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

bool EnthalpyTIF::distance1D_zold(const int i, const int j, const int k,
                                  const int dir, real & tint, real & dist) {
  real centrez = phi.zc(k);
  if(dir>0) {
    real edgez= centrez+0.5*phi.dzc(k);
    real intz = fsold[Comp::k()][i][j][k+1];
    real offz = phi.zc(k+1);
    if(intz>=centrez&&intz<=edgez) {
      tint = Tint_old(i,j,k);
      dist = intz-centrez;
      return true;
    } else if(intz>=edgez&&intz<=offz) {
      tint = Tint_old(i,j,k+1);
      dist = intz-centrez;
      return true;
    }
  } else {
    real edgez= centrez-0.5*phi.dzc(k);
    real intz = fsold[Comp::k()][i][j][k  ];
    real offz = phi.zc(k-1);
    if(intz<=centrez&&intz>=edgez) {
      tint = Tint_old(i,j,k);
      dist = centrez-intz;
      return true;
    } else if(intz>=offz&&intz<=edgez) {
      tint = Tint_old(i,j,k-1);
      dist = centrez-intz;
      return true;
    } 
  }
 
  return false;
}
