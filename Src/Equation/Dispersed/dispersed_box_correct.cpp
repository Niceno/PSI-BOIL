#include "dispersed.h"

/******************************************************************************/
void Dispersed::box_correct(int p) {

  /*=======================================================================+
  |                   store boxes in local coordinates.                    |
  |  If processor's domain is fully outside the box then put box indices   |
  |  to OFF.                                                               |  
  |  If processor's domain is fully inside the box the put box indices     |
  |  to processors's indices.                                              | 
  |  If processor's domain intersects with the box then box indices range  |
  |  up to processor indice.                                               |
  +=======================================================================*/
  const real xp = particles[p].x();
  const real yp = particles[p].y();
  const real zp = particles[p].z();
  const real rp = 0.5 * particles[p].d();

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

  /* store boxes in local coordinates */
  if(particles[p].si() != OFF || particles[p].ei() != OFF) {
    if(particles[p].si() == OFF) particles[p].si( si() );  
    if(particles[p].ei() == OFF) particles[p].ei( ei() );
  }
  if(particles[p].si() == OFF && particles[p].ei() == OFF) { 
    particles[p].si( si() );
    if(xc(si()) > xbw && xc(si()) < xbe) particles[p].ei( ei() ); 
  }
  if(particles[p].sj() != OFF || particles[p].ej() != OFF) {
    if(particles[p].sj() == OFF) particles[p].sj( sj() );  
    if(particles[p].ej() == OFF) particles[p].ej( ej() );
  }
  if(particles[p].sj() == OFF && particles[p].ej() == OFF) { 
    particles[p].sj( sj() );
    if(yc(sj()) > ybs && yc(sj()) < ybn) particles[p].ej( ej() ); 
  }
  if(particles[p].sk() != OFF || particles[p].ek() != OFF) {
    if(particles[p].sk() == OFF) particles[p].sk( sk() );  
    if(particles[p].ek() == OFF) particles[p].ek( ek() );
  }
  if(particles[p].sk() == OFF && particles[p].ek() == OFF) {
    particles[p].sk( sk() );
    if(zc(sk()) > zbb && zc(sk()) < zbt) particles[p].ek( ek() ); 
  }

}

/*----------------------------------------------------------------------------+
 '$Id: dispersed_box_correct.cpp,v 1.1 2015/08/19 11:29:18 badreddine Exp $'/
+-----------------------------------------------------------------------------*/
