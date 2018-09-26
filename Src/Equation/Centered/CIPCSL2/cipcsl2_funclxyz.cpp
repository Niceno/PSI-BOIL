#include "cipcsl2.h"
#include <cmath>

/******************************************************************************/
void CIPCSL2::funclxyz(real & outd
          ,real & outf
          ,real fi
          ,real fiup
          ,real d
          ,real den
          ,real x
          ,real isgn) {
#ifdef RCIP
    real s,bb,alpha,a,b,c;
    s=-isgn*den/d;
    real fiup_s = std::max(1e-12,fabs(fiup-s));
    bb=(fabs(s-fi)/fiup_s -1.0 )/d;
    if((s-fi)*(fiup-s)>-1.0e-16) {
      alpha=1.0;
    } else {
      alpha=0.0;
    }
    a=(fi-s+(fiup-s)*(1.0+alpha*bb*d))/(d*d);
    b=s*alpha*bb+(s-fi)/d-a*d;
    c=fi;

    outd=-(((a*x+b)*x+c)*x)/(1.0+alpha*bb*x);
    outf=((3.0*a*x+2.0*b)*x+c)/(1.0+alpha*bb*x)
        +outd*alpha*bb/(1.0+alpha*bb*x);
#else
    real a1,a2;
    a1=(fi+fiup)/(d*d)+2.0*isgn*den/(d*d*d);
    a2=-(2.0*fi+fiup)/d-3.0*isgn*den/(d*d);

    outf=3.0*a1*x*x+2.0*a2*x+fi;
    outd=-((a1*x+a2)*x+fi)*x ;
#endif
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_funclxyz.cpp,v 1.4 2011/11/07 08:17:10 sato Exp $'/
+-----------------------------------------------------------------------------*/
