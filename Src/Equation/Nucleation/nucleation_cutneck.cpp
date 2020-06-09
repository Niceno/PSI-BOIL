#include "nucleation.h"

/***************************************************************************//**
*  cut neck of the bubble when the bottom-area < area_neck
*******************************************************************************/
void Nucleation::cutneck (const real area_neck) {

  /* set neck */
  for(int ns=0; ns<size(); ns++){
    bool bneck=false;
    real a_vapor=0.0;
    if( time->current_time() > (sites[ns].time_seed() + seed_period) ) {
      a_vapor
	= area_vapor_sum( Range<real>(sites[ns].x()-rcut, sites[ns].x()+rcut)
                        , Range<real>(sites[ns].y()-rcut, sites[ns].y()+rcut)
                        , Range<real>(zbtm, zbtm+dxmin));

      /*--------------------------------------------------------+
      |  cut neck if (1) vapor area is smaller than area_neck   |
      |  and (2) replant process is not running                 |
      +------------------------------------------------------- */
      if( a_vapor < area_neck &&
          !sites[ns].seed() ) {
        bneck=true;

        /* store the beginning time of cutneck */
        if( !sites[ns].neck_prev() ) {
          sites[ns].set_time_cutneck( time->current_time() );
        }

      }

    }
    sites[ns].set_neck(bneck);
    sites[ns].set_neck_prev(bneck);
    boil::oout<<"cutneck: "<<time->current_time()<<" "<<ns<<" "<<bneck<<" "
              <<" "<<a_vapor<<" "<<sites[ns].neck_prev()<<"\n";
  }
  /* copy bneck to dummy sites */
  for(int nsd=0; nsd<dsize(); nsd++){
    int ns=dsites[nsd].father();
    bool bneck = sites[ns].neck();
    dsites[nsd].set_neck(bneck);
  }

  /* cut neck */
  for(int ns=0; ns<size(); ns++){
    if( sites[ns].neck() ){
      real xfirst = sites[ns].x() - rcut;
      real xlast  = sites[ns].x() + rcut;
      real yfirst = sites[ns].y() - rcut;
      real ylast  = sites[ns].y() + rcut;
      real zfirst = zbtm;
      real zlast  = zbtm+dxmin;
      for (int i=(*clr).si(); i<=(*clr).ei(); i++) {
        if ((*clr).xc(i)<xfirst) continue;
        if ((*clr).xc(i)>xlast ) continue;
        for (int j=(*clr).sj(); j<=(*clr).ej(); j++) {
          if ((*clr).yc(j)<yfirst) continue;
          if ((*clr).yc(j)>ylast ) continue;
	  for (int k=(*clr).sk(); k<=(*clr).ek(); k++) {
            if ((*clr).zc(k)<zfirst) continue;
            if ((*clr).zc(k)>zlast ) continue;
            (*vf)[i][j][k]= (sig>0) ? 1.0 : 0.0;
	    //std::cout<<i<<" "<<j<<" "<<k<<"\n";
          }
        }
      }
    }
  }

  /* dummy sites */
  for(int nsd=0; nsd<dsize(); nsd++){
    if( dsites[nsd].neck() ){
      real xfirst = dsites[nsd].x() - rcut;
      real xlast  = dsites[nsd].x() + rcut;
      real yfirst = dsites[nsd].y() - rcut;
      real ylast  = dsites[nsd].y() + rcut;
      real zfirst = zbtm;
      real zlast  = zbtm+dxmin;
      for (int i=(*clr).si(); i<=(*clr).ei(); i++) {
        if ((*clr).xc(i)<xfirst) continue;
        if ((*clr).xc(i)>xlast ) continue;
        for (int j=(*clr).sj(); j<=(*clr).ej(); j++) {
          if ((*clr).yc(j)<yfirst) continue;
          if ((*clr).yc(j)>ylast ) continue;
          for (int k=(*clr).sk(); k<=(*clr).ek(); k++) {
            if ((*clr).zc(k)<zfirst) continue;
            if ((*clr).zc(k)>zlast ) continue;
            (*vf)[i][j][k]= (sig>0) ? 1.0 : 0.0;
          }
        }
      }
    }
  }
}
