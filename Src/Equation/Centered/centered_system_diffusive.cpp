#include "centered.h"

/***************************************************************************//**
*  \brief Creates diffusive part of the system matrix \f$ [A] \f$.
*******************************************************************************/
void Centered::create_system_diffusive(const Property * f_prop,   
                                       const Property * s_prop,
                                       const Scalar * diff_eddy) {

  /* initialize: get time stepping coefficient */
  real tsc = diff_ts.N();
  assert( tsc > 0.0 );

  /*------------------------------------+
  |                                     |
  |  no conduction through solid parts  |
  |                                     |
  +------------------------------------*/
  if( !solid() ) {

    /* coefficients in i direction (w and e) */
    for_ijk(i,j,k) {
      const real a_x = dSx(i,j,k);
      real lm = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i-1,j,k) );
      real lp = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i+1,j,k) );
      if (diff_eddy) {
        lm += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i-1][j][k]);
        lp += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i+1][j][k]);
      }
      A.w[i][j][k] = tsc * lm * a_x / dxw(i);
      A.e[i][j][k] = tsc * lp * a_x / dxe(i);
    }

    /* coefficients in j direction (s and n) */
    for_ijk(i,j,k) {
      const real a_y = dSy(i,j,k);
      real lm = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j-1,k) );
      real lp = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j+1,k) );
      if (diff_eddy) {
        lm += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j-1][k]);
        lp += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j+1][k]);
      }
      A.s[i][j][k] = tsc * lm * a_y / dys(j);
      A.n[i][j][k] = tsc * lp * a_y / dyn(j);
    }
  
    /* coefficients in k direction (b and t) */
    for_ijk(i,j,k) {
      const real a_z = dSz(i,j,k);
      real lm = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j,k-1) );
      real lp = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j,k+1) );
      if (diff_eddy) {
        lm += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j][k-1]);
        lp += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j][k+1]);
      }
      A.b[i][j][k] = tsc * lm * a_z / dzb(k);
      A.t[i][j][k] = tsc * lp * a_z / dzt(k);
    }

    /*-------------------------------+
    |  a "touch" from immersed body  |
    +-------------------------------*/
    if(dom->ibody().nccells() > 0) {
      for(int cc=0; cc<dom->ibody().nccells(); cc++) {
        int i,j,k;
        dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);
  
        /* w */
        if( dom->ibody().on(i-1,j,k) ) {
          const real fSw  = dom->ibody().fSw(cc);
          const real fdxw = dom->ibody().fdxw(cc);
          if( dom->ibody().on(i,j,k) ) A.w[i]  [j][k] *= (fSw / fdxw);
          else if(fdxw != 1.0)         A.e[i-1][j][k] /= (1.0 - fdxw);
        }
        
        /* e */
        if( dom->ibody().on(i+1,j,k) ) {
          const real fSe  = dom->ibody().fSe(cc);
          const real fdxe = dom->ibody().fdxe(cc);
          if( dom->ibody().on(i,j,k) ) A.e[i]  [j][k] *= (fSe / fdxe);
          else if(fdxe != 1.0)         A.w[i+1][j][k] /= (1.0 - fdxe);
        }
  
        /* s */
        if( dom->ibody().on(i,j-1,k) ) {
          const real fSs  = dom->ibody().fSs(cc);
          const real fdys = dom->ibody().fdys(cc);
          if( dom->ibody().on(i,j,k) ) A.s[i][j]  [k] *= (fSs / fdys);
          else if(fdys != 1.0)         A.n[i][j-1][k] /= (1.0 - fdys);
        }
          
        /* n */
        if( dom->ibody().on(i,j+1,k) ) {
          const real fSn  = dom->ibody().fSn(cc);
          const real fdyn = dom->ibody().fdyn(cc);
          if( dom->ibody().on(i,j,k) ) A.n[i][j]  [k] *= (fSn / fdyn);
          else if(fdyn != 1.0)         A.s[i][j+1][k] /= (1.0 - fdyn);
        }
  
        /* b */
        if( dom->ibody().on(i,j,k-1) ) {
          const real fSb  = dom->ibody().fSb(cc);
          const real fdzb = dom->ibody().fdzb(cc);
          if( dom->ibody().on(i,j,k) ) A.b[i][j][k]   *= (fSb / fdzb);
          else if(fdzb != 1.0)         A.t[i][j][k-1] /= (1.0 - fdzb);
        }
  
        /* t */
        if( dom->ibody().on(i,j,k+1) ) {
          const real fSt  = dom->ibody().fSt(cc);
          const real fdzt = dom->ibody().fdzt(cc);
          if( dom->ibody().on(i,j,k) ) A.t[i][j][k]   *= (fSt / fdzt);
          else if(fdzt != 1.0)         A.b[i][j][k+1] /= (1.0 - fdzt);
        }

      }
  
    } /* is there an immersed body */

  /*------------------------------------------+
  |                                           |
  |  features conduction through solid parts  |
  |                                           |
  +------------------------------------------*/
  } else {

    /* coefficients in i direction (w and e) */
    for_ijk(i,j,k) {
      real a_x = dSx(i,j,k);
      real a_w = dSx(i,j,k);
      real a_e = dSx(i,j,k);
      if(dom->ibody().cut(i,j,k)) {
        a_w *= dom->ibody().fSw(i,j,k);
        a_e *= dom->ibody().fSe(i,j,k);
      }
      real lmf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i-1,j,k) );
      real lpf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i+1,j,k) );
      if (diff_eddy) {
        lmf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i-1][j][k]);
        lpf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i+1][j][k]);
      }
      const real lms = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i-1,j,k));
      const real lps = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i+1,j,k));
      /* w */
      if( dom->ibody().on (i-1,j,k) && dom->ibody().on (i,j,k) )  
        A.w[i][j][k] = tsc * ( lmf*a_w + lms*(a_x-a_w) ) / dxw(i);
      else if( dom->ibody().off(i-1,j,k) && dom->ibody().off(i,j,k) )  
        A.w[i][j][k] = tsc * ( lms*a_w + lmf*(a_x-a_w) ) / dxw(i);
      else {
        real lamc, lamn, dc, dn, xw_surf;
        if( dom->ibody().fdxw(i,j,k) < 1.0 )
          xw_surf = phi.xc(i)   - dom->ibody().fdxw(i,  j,k) * dxw(i);
        else
          xw_surf = phi.xc(i-1) + dom->ibody().fdxe(i-1,j,k) * dxe(i-1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = phi.xc(i) - xw_surf;
        dn = xw_surf - phi.xc(i-1);
        A.w[i][j][k] = tsc * a_x * lamc*lamn / (lamc*dn + lamn*dc);
      }
      
      /* e */
      if( dom->ibody().on (i+1,j,k) && dom->ibody().on (i,j,k) )  
        A.e[i][j][k] = tsc * ( lpf*a_e + lps*(a_x-a_e) ) / dxe(i);
      else if( dom->ibody().off(i+1,j,k) && dom->ibody().off(i,j,k) )  
        A.e[i][j][k] = tsc * ( lps*a_e + lpf*(a_x-a_e) ) / dxe(i);
      else {
        real lamc, lamn, dc, dn, xe_surf;
        if( dom->ibody().fdxe(i,j,k) < 1.0 )
          xe_surf = phi.xc(i)   + dom->ibody().fdxe(i,  j,k) * dxe(i);
        else
          xe_surf = phi.xc(i+1) - dom->ibody().fdxw(i+1,j,k) * dxw(i+1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = xe_surf - phi.xc(i);
        dn = phi.xc(i+1) - xe_surf;
        A.e[i][j][k] = tsc * a_x * lamc*lamn / (lamc*dn + lamn*dc);
      }
    }

    /* coefficients in j direction (s and n) */
    for_ijk(i,j,k) {
      real a_y = dSy(i,j,k);
      real a_s = dSy(i,j,k);
      real a_n = dSy(i,j,k);
      if(dom->ibody().cut(i,j,k)) {
        a_s *= dom->ibody().fSs(i,j,k);
        a_n *= dom->ibody().fSn(i,j,k);
      }
      real lmf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j-1,k) );
      real lpf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j+1,k) );
      if (diff_eddy) {
        lmf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j-1][k]);
        lpf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j+1][k]);
      }
      const real lms = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i,j-1,k));
      const real lps = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i,j+1,k));
      /* s */
      if( dom->ibody().on (i,j-1,k) && dom->ibody().on (i,j,k) )  
        A.s[i][j][k] = tsc * ( lmf*a_s + lms*(a_y-a_s) ) / dys(j);
      else if( dom->ibody().off(i,j-1,k) && dom->ibody().off(i,j,k) ) 
        A.s[i][j][k] = tsc * ( lms*a_s + lmf*(a_y-a_s) ) / dys(j);
      else {
        real lamc, lamn, dc, dn, ys_surf;
        if( dom->ibody().fdys(i,j,k) < 1.0 )
          ys_surf = phi.yc(j)   - dom->ibody().fdys(i,j,  k) * dys(j);
        else
          ys_surf = phi.yc(j-1) + dom->ibody().fdyn(i,j-1,k) * dyn(j-1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = phi.yc(j) - ys_surf;
        dn = ys_surf - phi.yc(j-1);
        A.s[i][j][k] = tsc * a_y * lamc*lamn / (lamc*dn + lamn*dc);
      }
      
      /* n */
      if( dom->ibody().on (i,j+1,k) && dom->ibody().on (i,j,k) )  
        A.n[i][j][k] = tsc * ( lpf*a_n + lps*(a_y-a_n) ) / dyn(j);
      else if( dom->ibody().off(i,j+1,k) && dom->ibody().off(i,j,k) )  
        A.n[i][j][k] = tsc * ( lps*a_n + lpf*(a_y-a_n) ) / dyn(j);
      else {
        real lamc, lamn, dc, dn, yn_surf;
        if( dom->ibody().fdyn(i,j,k) < 1.0 )
          yn_surf = phi.yc(j)   + dom->ibody().fdyn(i,j,  k) * dyn(j);
        else
          yn_surf = phi.yc(j+1) - dom->ibody().fdys(i,j+1,k) * dys(j+1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = yn_surf - phi.yc(j);
        dn = phi.yc(j+1) - yn_surf;
        A.n[i][j][k] = tsc * a_y * lamc*lamn / (lamc*dn + lamn*dc);
      }
    }
  
    /* coefficients in k direction (b and t) */
    for_ijk(i,j,k) {
      real a_z = dSz(i,j,k);
      real a_b = dSz(i,j,k);
      real a_t = dSz(i,j,k);
      if(dom->ibody().cut(i,j,k)) {
        a_b *= dom->ibody().fSb(i,j,k);
        a_t *= dom->ibody().fSt(i,j,k);
      }
      real lmf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j,k-1) );
      real lpf = 0.5 * ( f_prop->value(i,j,k) + f_prop->value(i,j,k+1) );
      if (diff_eddy) {
        lmf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j][k-1]);
        lpf += 0.5 * ((*diff_eddy)[i][j][k] + (*diff_eddy)[i][j][k+1]);
      }
      const real lms = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i,j,k-1));
      const real lps = 0.5 * (solid()->lambda(i,j,k)+solid()->lambda(i,j,k+1));
      /* b */
      if( dom->ibody().on(i,j,k-1) && dom->ibody().on(i,j,k) ) 
        A.b[i][j][k] = tsc * ( lmf*a_b + lms*(a_z-a_b) ) / dzb(k);
      else if( dom->ibody().off(i,j,k-1) && dom->ibody().off(i,j,k) ) 
        A.b[i][j][k] = tsc * ( lms*a_b + lmf*(a_z-a_b) ) / dzb(k);
      else {
        real lamc, lamn, dc, dn, zb_surf;
        if( dom->ibody().fdzb(i,j,k) < 1.0 )
          zb_surf = phi.zc(k)   - dom->ibody().fdzb(i,j,k)   * dzb(k);
        else
          zb_surf = phi.zc(k-1) + dom->ibody().fdzt(i,j,k-1) * dzt(k-1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = phi.zc(k) - zb_surf;
        dn = zb_surf - phi.zc(k-1);
        A.b[i][j][k] = tsc * a_z * lamc*lamn / (lamc*dn + lamn*dc);
      }
      /* t */
      if( dom->ibody().on (i,j,k+1) && dom->ibody().on (i,j,k) ) 
        A.t[i][j][k] = tsc * ( lpf*a_t + lps*(a_z-a_t) ) / dzt(k);
      else if( dom->ibody().off(i,j,k+1) && dom->ibody().off(i,j,k) ) 
        A.t[i][j][k] = tsc * ( lps*a_t + lpf*(a_z-a_t) ) / dzt(k);
      else {
        real lamc, lamn, dc, dn, zt_surf;
        if( dom->ibody().fdzt(i,j,k) < 1.0 )
          zt_surf = phi.zc(k)   + dom->ibody().fdzt(i,j,k)   * dzt(k);
        else
          zt_surf = phi.zc(k+1) - dom->ibody().fdzb(i,j,k+1) * dzb(k+1);
        if( dom->ibody().on(i,j,k) ) {
          lamc = f_prop->value(i,j,k);
          if (diff_eddy) lamc += (*diff_eddy)[i][j][k];
          lamn = solid()->lambda(i,j,k);
        } else {
          lamc = solid()->lambda(i,j,k);
          lamn = f_prop->value(i,j,k);
          if (diff_eddy) lamn += (*diff_eddy)[i][j][k];
        }
        dc = zt_surf - phi.zc(k);
        dn = phi.zc(k+1) - zt_surf;
        A.t[i][j][k] = tsc * a_z * lamc*lamn / (lamc*dn + lamn*dc);
      }
    }

  } /* conduction through solid */
 
  /*---------------------------------+
  |  summ the central coefficent up  |
  +---------------------------------*/
  for_ijk(i,j,k) {
    A.c[i][j][k] += A.w[i][j][k] + A.e[i][j][k]
                  + A.s[i][j][k] + A.n[i][j][k]
                  + A.b[i][j][k] + A.t[i][j][k];
  }

  /*---------------------+ 
  |  this is needed too  |
  +---------------------*/
  if( !solid() ) 
    if(dom->ibody().nccells() > 0) {
      for_ijk(i,j,k)
        if( dom->ibody().off(i,j,k) ) {
          A.c[i][j][k]  = 1.0;          
          A.w[i][j][k]  = 0.0;
          A.e[i][j][k]  = 0.0;
          A.s[i][j][k]  = 0.0;
          A.n[i][j][k]  = 0.0;
          A.b[i][j][k]  = 0.0;
          A.t[i][j][k]  = 0.0;
          A.ci[i][j][k] = 1.0;
          fnew[i][j][k] = 0.0;
        }
    } /* is there an immersed body */

  A.c.exchange();
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();
  A.ci.exchange();

}
