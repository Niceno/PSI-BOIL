#include "scalar.h"

/******************************************************************************/
/*  interpolation using inverse of distance                                   */
/*  (x, y, z) does not have to locate in the decomposed domain.               */
/*  If it is outside of domain, it returns -boil::unreal.                     */
/*                                                                            */
/*  WARNING:: everytime mpi data exchange happens.                            */
/******************************************************************************/
real Scalar::Interpolate(const real x, const real y, const real z) const {
  if (x < dom->global_min_x() || x > dom->global_max_x() ||
      y < dom->global_min_y() || y > dom->global_max_y() ||
      z < dom->global_min_z() || z > dom->global_max_z() ) {
     boil::oout<<"Scalar_Interpolate: WARNING!!! ("<<x<<", "<<y<<", "<<z
                <<") is outside of computational domain.\n";
     return -boil::unreal;
  } else {
    real val_interpolate = interpolate(x, y, z);
    boil::cart.max_real(&val_interpolate);
    return val_interpolate;
  }
}
/******************************************************************************/
/*  3D linear interpolation                                                   */
/*  (x, y, z) must locate in the decomposed domain, otherwise it returns      */
/*  -boil::unreal                                                             */
/******************************************************************************/
real Scalar::interpolate(const real x, const real y, const real z) const {

  if (x < dom->global_min_x() || x > dom->global_max_x() ||
      y < dom->global_min_y() || y > dom->global_max_y() ||
      z < dom->global_min_z() || z > dom->global_max_z() ) {
     boil::oout<<"Scalar_Interpolate: WARNING!!! ("<<x<<", "<<y<<", "<<z
                <<") is outside of computational domain.\n";
     return -boil::unreal;
  }

  real val_interpolate = -boil::unreal;  // default return value
					 // if (x,y,z) is not included in the
					 // decomposed domain, this value will
					 // be returned.

  if( dom->contains_xyz(x, y, z) ){ // (x, y, z) is inside of decomposed domain
    int ic = i(x);  // cell index i, which includes x
    int jc = j(y);  // cell index j, which includes y
    int kc = k(z);  // cell index k, which includes z
    if (ic<si() || ic>ei() || jc<sj() || jc>ej() || kc<sk() || kc>ek()) {
      std::cout<<"Error:scalar_interpolate:proc= "<<boil::cart.iam()
	       <<" ic,jc,kc="<<ic<<" "<<jc<<" "<<kc
               <<" si,ei "<<si()<<" "<<ei()
               <<" sj,ej "<<sj()<<" "<<ej()
               <<" sk,ek "<<sk()<<" "<<ek()<<"\n";
      exit(0);
    }
    // The point (x, y, z) locates within the region in the cell centers
    // [im][jm][km] and [ip][jp][kp].
    // The value at (x,y,z) is interpolated from the eight points.
    int im = ic - 1;
    int jm = jc - 1;
    int km = kc - 1;
    if (x > xc(ic) ) im = ic;
    if (y > yc(jc) ) jm = jc;
    if (z > zc(kc) ) km = kc;

    int ip = im +1;
    int jp = jm +1;
    int kp = km +1;

#if 0
    // inverse of distance
    real w[2][2][2];  // weight function
    real sum_w = 0.0;
    real eps = 1.0e-8 * std::min(xc(ic), std::min(yc(jc), zc(kc)));
    for (int ii = im; ii <= ip; ii++) {
      for (int jj = jm; jj <= jp; jj++) {
        for (int kk = km; kk <= kp; kk++) {
          real dist = sqrt( pow(x-xc(ii),2) + pow(y-yc(jj),2) + pow(z-zc(kk),2));
          dist = std::max(dist, eps);
          real disti = 1.0/dist;
          w[ii-im][jj-jm][kk-km] = disti;
	  sum_w += disti;
        }
      }
    }
    val_interpolate = 0.0;
    for (int ii = im; ii <= ip; ii++) {
      for (int jj = jm; jj <= jp; jj++) {
        for (int kk = km; kk <= kp; kk++) {
           val_interpolate += w[ii-im][jj-jm][kk-km] / sum_w * val[ii][jj][kk];
        } 
      }
    }
#endif
#if 1
  // 3D linear interpolation
  real x_nrm = (x-xc(im))/(xc(ip) -xc(im));  // normalized x
  real y_nrm = (y-yc(jm))/(yc(jp) -yc(jm));  // normalized y
  real z_nrm = (z-zc(km))/(zc(kp) -zc(km));  // normalized z
  val_interpolate = val[im][jm][km]*(1.0-x_nrm)*(1.0-y_nrm)*(1.0-z_nrm)
                  + val[ip][jm][km]*(    x_nrm)*(1.0-y_nrm)*(1.0-z_nrm)
                  + val[im][jp][km]*(1.0-x_nrm)*(    y_nrm)*(1.0-z_nrm)
                  + val[ip][jp][km]*(    x_nrm)*(    y_nrm)*(1.0-z_nrm)
                  + val[im][jm][kp]*(1.0-x_nrm)*(1.0-y_nrm)*(    z_nrm)
                  + val[ip][jm][kp]*(    x_nrm)*(1.0-y_nrm)*(    z_nrm) 
                  + val[im][jp][kp]*(1.0-x_nrm)*(    y_nrm)*(    z_nrm)
                  + val[ip][jp][kp]*(    x_nrm)*(    y_nrm)*(    z_nrm);
#endif
  }

  return val_interpolate;
}

