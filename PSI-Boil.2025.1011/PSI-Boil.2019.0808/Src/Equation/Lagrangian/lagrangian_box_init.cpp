#include "lagrangian.h"

/******************************************************************************/
void Lagrangian::box_init(int p) {

  /*=======================================================+
  |                                                        |
  |  browse through all particles to store bounding boxes  |
  |                                                        |
  +=======================================================*/

  /* create a box around the particle */
  const real xp = particles[p].x();
  const real yp = particles[p].y();
  const real zp = particles[p].z();
  const real rp = 0.5 * particles[p].d();

  /* check if particle is inside the computational domain */
  /*assert( xp >= dom->global_min_x() );
  assert( xp <= dom->global_max_x() );
  assert( yp >= dom->global_min_y() );
  assert( yp <= dom->global_max_y() );
  assert( zp >= dom->global_min_z() );
  assert( zp <= dom->global_max_z() );*/

int ind = 0;
if( (xp >= dom->global_min_x()) && (xp <= dom->global_max_x()) ) {ind=1;} else{ind=0;};
if( (yp >= dom->global_min_y()) && (yp <= dom->global_max_y()) ) {ind=1;} else{ind=0;};
if( (zp >= dom->global_min_z()) && (zp <= dom->global_max_z()) ) {ind=1;} else{ind=0;};

/*if( (xp >= dom->global_min_x()) &&
    (xp <= dom->global_max_x()) &&
    (yp >= dom->global_min_y()) &&
    (yp <= dom->global_max_y()) &&
    (zp >= dom->global_min_z()) &&
    (zp <= dom->global_max_z()) ) {*/  //can lead to error-Segmentation fault (core dumped)

if(ind=1) {
  /* create bounding box coordinates (w, e, s, n, b, t),
     taking care of the extents of the computational domain too */

  const real dx = dxc(2);
  const real dy = dyc(2); 
  const real dz = dzc(2);
  const real bdr = box_diam_ratio;

  const real xbw = boil::maxr(xp-rp * bdr - dx,
                   dom->global_min_x() + boil::pico);
  const real xbe = boil::minr(xp+rp * bdr + dx,
                   dom->global_max_x() - boil::pico);
  const real ybs = boil::maxr(yp-rp * bdr - dy,
                   dom->global_min_y() + boil::pico);
  const real ybn = boil::minr(yp+rp * bdr + dy,
                   dom->global_max_y() - boil::pico);
  const real zbb = boil::maxr(zp-rp * bdr - dz,
                   dom->global_min_z() + boil::pico);
  const real zbt = boil::minr(zp+rp * bdr + dz,
                   dom->global_max_z() - boil::pico);

  /* store box dimensions (needed for gradients later on) */
  particles[p].box_dx(xbe - xbw);
  particles[p].box_dy(ybn - ybs);
  particles[p].box_dz(zbt - zbb);

  /* store boxes in local coordinates */
  particles[p].si( dom->local_i( dom->I(xbw) ) );  
  particles[p].ei( dom->local_i( dom->I(xbe) ) );  
  particles[p].sj( dom->local_j( dom->J(ybs) ) );  
  particles[p].ej( dom->local_j( dom->J(ybn) ) ); 
  particles[p].sk( dom->local_k( dom->K(zbb) ) );  
  particles[p].ek( dom->local_k( dom->K(zbt) ) ); 
  }

  //else erase(p);

}

/*----------------------------------------------------------------------------+
 '$Id: lagrangian_box_init.cpp,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
