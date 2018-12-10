#include "finescalar.h"

/* calculate the cross product of the arguments */
FineScalar::XYZ FineScalar::CrossProduct(const XYZ p1, const XYZ p2) {
  XYZ cp;
  cp.x = p1.y * p2.z - p1.z * p2.y;
  cp.y = p1.z * p2.x - p1.x * p2.z;
  cp.z = p1.x * p2.y - p1.y * p2.x;

  return(cp);
}

/* calculate the dot product of the arguments */
real FineScalar::DotProduct(const XYZ p1, const XYZ p2) {
  real dp(0.0);
  dp += p1.x * p2.x;
  dp += p1.y * p2.y;
  dp += p1.z * p2.z;

  return(dp);
}

/* add two xyz vectors */
FineScalar::XYZ FineScalar::PlusXYZ(const XYZ p1, const XYZ p2) {
  XYZ pv;
  pv.x = p1.x + p2.x;
  pv.y = p1.y + p2.y;
  pv.z = p1.z + p2.z;

  return(pv);
}

real FineScalar::calc_area(const std::vector<XYZ> &vect, const XYZ norm) { 
  XYZ cross;
  cross.x = 0.0;
  cross.y = 0.0;
  cross.z = 0.0;
  for(int i = 0; i != vect.size(); ++i) {
    cross = PlusXYZ(cross,CrossProduct(vect[i],vect[(i+1) % vect.size()]));
  }
  return 0.5 * fabs(DotProduct(cross,norm));
}

/******************************************************************************/
void FineScalar::cal_adens_geom() {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adens, adensgeom
*******************************************************************************/

  boil::timer.start("finescalar adens geom");

  /* cell centered */
  for_vijk(adens,i,j,k) {
    real gradient = 0.0;
   
    real dx = adens.dxc(i);
    real dy = adens.dyc(j);
    real dz = adens.dzc(k);
 
    /* x-gradient */
    real marker_w   = value(i,j,k, w  ()); 
    real marker_ws  = value(i,j,k, ws ()); 
    real marker_wn  = value(i,j,k, wn ()); 
    real marker_wb  = value(i,j,k, wb ()); 
    real marker_wt  = value(i,j,k, wt ()); 
    real marker_wsb = value(i,j,k, wsb()); 
    real marker_wst = value(i,j,k, wst()); 
    real marker_wnb = value(i,j,k, wnb()); 
    real marker_wnt = value(i,j,k, wnt()); 

    real marker_e   = value(i,j,k, e  ()); 
    real marker_es  = value(i,j,k, es ()); 
    real marker_en  = value(i,j,k, en ()); 
    real marker_eb  = value(i,j,k, eb ()); 
    real marker_et  = value(i,j,k, et ()); 
    real marker_esb = value(i,j,k, esb()); 
    real marker_est = value(i,j,k, est()); 
    real marker_enb = value(i,j,k, enb()); 
    real marker_ent = value(i,j,k, ent()); 

    real gradx = grad_1D(marker_w,
                         marker_ws, marker_wn, marker_wb, marker_wt,
                         marker_wsb, marker_wst, marker_wnb, marker_wnt,
                         marker_e,
                         marker_es, marker_en, marker_eb, marker_et,
                         marker_esb, marker_est, marker_enb, marker_ent,
                         dx, k);

    /* y-gradient */
    real marker_s   = value(i,j,k, s  ()); 
    real marker_sw  = marker_ws; 
    real marker_se  = marker_es; 
    real marker_sb  = value(i,j,k, sb ()); 
    real marker_st  = value(i,j,k, st ()); 
    real marker_swb = marker_wsb; 
    real marker_swt = marker_wst; 
    real marker_seb = marker_esb; 
    real marker_set = marker_est; 

    real marker_n   = value(i,j,k, n  ()); 
    real marker_nw  = marker_wn; 
    real marker_ne  = marker_en; 
    real marker_nb  = value(i,j,k, nb ()); 
    real marker_nt  = value(i,j,k, nt ()); 
    real marker_nwb = marker_wnb; 
    real marker_nwt = marker_wnt; 
    real marker_neb = marker_enb; 
    real marker_net = marker_ent; 

    real grady = grad_1D(marker_s,
                         marker_sw, marker_se, marker_sb, marker_st,
                         marker_swb, marker_swt, marker_seb, marker_set,
                         marker_n,
                         marker_nw, marker_ne, marker_nb, marker_nt,
                         marker_nwb, marker_nwt, marker_neb, marker_net,
                         dy, k);

    /* z-gradient */
    real marker_b   = value(i,j,k, b  ()); 
    real marker_bw  = marker_wb; 
    real marker_be  = marker_eb; 
    real marker_bs  = marker_sb; 
    real marker_bn  = marker_nb; 
    real marker_bws = marker_wsb; 
    real marker_bwn = marker_wnb; 
    real marker_bes = marker_esb; 
    real marker_ben = marker_enb; 

    real marker_t   = value(i,j,k, t  ()); 
    real marker_tw  = marker_wt; 
    real marker_te  = marker_et; 
    real marker_ts  = marker_st; 
    real marker_tn  = marker_nt; 
    real marker_tws = marker_wst; 
    real marker_twn = marker_wnt; 
    real marker_tes = marker_est; 
    real marker_ten = marker_ent; 

    real gradz = grad_1D(marker_b,
                         marker_bw, marker_be, marker_bs, marker_bn, 
                         marker_bws, marker_bwn, marker_bes, marker_ben,
                         marker_t,
                         marker_tw, marker_te, marker_ts, marker_tn,
                         marker_tws, marker_twn, marker_tes, marker_ten,
                         dz, k);

    adens[i][j][k] = sqrt(gradx*gradx+grady*grady+gradz*gradz);

    if(adens[i][j][k]>boil::pico) {
      
      /* calculate normal vector at cell center */
#if 0
      /* m points to the gas, normalized */
      real m1 = -gradx;
      real m2 = -grady;
      real m3 = -gradz;
       
      real mmag = sqrt(m1*m1+m2*m2+m3*m3);

      m1 /= mmag;
      m2 /= mmag;
      m3 /= mmag;

      /* vn points to the gas, normalized, in normalized space */
      real vn1 = m1*dx;
      real vn2 = m2*dy;
      real vn3 = m3*dz;

      real vnmag = sqrt(vn1*vn1+vn2*vn2+vn3*vn3);

      vn1 /= vnmag;
      vn2 /= vnmag;
      vn3 /= vnmag;
#else
      /* vn points to the gas, normalized, in normalized space */
      real vn1 = -(*nx)[i][j][k];
      real vn2 = -(*ny)[i][j][k];
      real vn3 = -(*nz)[i][j][k];
      /* m points to the gas, normalized */
      real m1 = vn1/dx;
      real m2 = vn2/dy;
      real m3 = vn3/dz;

      real mmag = sqrt(m1*m1+m2*m2+m3*m3);
      m1 /= mmag;
      m2 /= mmag;
      m3 /= mmag;

      real vnmag = 1.0/mmag;
#endif

      /* vm is normalized and standardized */
      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      real vm3 = fabs(vn3);

      real denom = vm1+vm2+vm3;
      real qa = 1.0/denom;

      /* alpha is standardized */
      real c = (*phi)[i][j][k]; 
      real alpha = denom*calc_alpha(c, vm1*qa, vm2*qa, vm3*qa);

      /* remove the effect of mirroring */
      if(m1<0.) alpha -= vm1;
      if(m2<0.) alpha -= vm2;
      if(m3<0.) alpha -= vm3;

      /* remove the effect of normalization */
      alpha *= vnmag;

#if 0
      /* translate alpha */
      real xc = (*phi).xc(i);
      real yc = (*phi).yc(j);
      real zc = (*phi).zc(k);

      alpha += m1*xc+m2*yc+m3*zc;
#endif

      XYZ normvector;
      normvector.x = m1;
      normvector.y = m2;
      normvector.z = m3;

      /* for each edge of the cell, try to find an intersection */
      std::vector<XYZ> pointset;

      bool xnorm_nonzero(vm1>boil::pico);
      bool ynorm_nonzero(vm2>boil::pico);
      bool znorm_nonzero(vm3>boil::pico);

      real xval, yval, zval;

      /* y = 0, z = 0, x varies */
      yval = 0.0;
      zval = 0.0;
 
      if(xnorm_nonzero) { 
        xval = (alpha-m2*yval-m3*zval)/m1;
        if(xval>=0.0&&xval<=dx) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval; 
          newpoint.z = zval;
          pointset.push_back(newpoint);
        } 
      }

      /* y = 0, z = dz, x varies */
      yval = 0.0;
      zval = dz;

      if(xnorm_nonzero) {
        xval = (alpha-m2*yval-m3*zval)/m1;
        if(xval>=0.0&&xval<=dx) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* y = dy, z = 0, x varies */
      yval = dy;
      zval = 0.0;

      if(xnorm_nonzero) {
        xval = (alpha-m2*yval-m3*zval)/m1;
        if(xval>=0.0&&xval<=dx) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* y = dy, z = dz, x varies */
      yval = dy;
      zval = dz;

      if(xnorm_nonzero) {
        xval = (alpha-m2*yval-m3*zval)/m1;
        if(xval>=0.0&&xval<=dx) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, z = 0, y varies */
      xval = 0.0;
      zval = 0.0;

      if(ynorm_nonzero) {
        yval = (alpha-m1*xval-m3*zval)/m2;
        if(yval>=0.0&&yval<=dy) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, z = dz, y varies */
      xval = 0.0;
      zval = dz;

      if(ynorm_nonzero) {
        yval = (alpha-m1*xval-m3*zval)/m2;
        if(yval>=0.0&&yval<=dy) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, z = 0, y varies */
      xval = dx;
      zval = 0.0;

      if(ynorm_nonzero) {
        yval = (alpha-m1*xval-m3*zval)/m2;
        if(yval>=0.0&&yval<=dy) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, z = dz, y varies */
      xval = dx;
      zval = dz;

      if(ynorm_nonzero) {
        yval = (alpha-m1*xval-m3*zval)/m2;
        if(yval>=0.0&&yval<=dy) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, y = 0, z varies */
      xval = 0.0;
      yval = 0.0;

      if(znorm_nonzero) {
        zval = (alpha-m1*xval-m2*yval)/m3;
        if(zval>=0.0&&zval<=dz) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, y = 0, z varies */
      xval = dx;
      yval = 0.0;

      if(znorm_nonzero) {
        zval = (alpha-m1*xval-m2*yval)/m3;
        if(zval>=0.0&&zval<=dz) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, y = dy, z varies */
      xval = 0.0;
      yval = dy;

      if(znorm_nonzero) {
        zval = (alpha-m1*xval-m2*yval)/m3;
        if(zval>=0.0&&zval<=dz) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }


      /* x = dx, y = dy, z varies */
      xval = dx;
      yval = dy;

      if(znorm_nonzero) {
        zval = (alpha-m1*xval-m2*yval)/m3;
        if(zval>=0.0&&zval<=dz) {
          XYZ newpoint;
          newpoint.x = xval;
          newpoint.y = yval;
          newpoint.z = zval;
          pointset.push_back(newpoint);
        }
      }

      if(pointset.size()<3) {
        adensgeom[i][j][k] = 0.0;
        //boil::oout<<"FS::caladensgeom: Warning, inconsistent geometry"<<boil::endl;
      } else {
        /* order points */
        /* for a convex polygon, neighboring vertices are the closest ones */
        std::vector<XYZ> orderedset;
        int numpoints = pointset.size();

        /* first point is arbitrary */  
        orderedset.push_back(pointset[0]);
        pointset.erase(pointset.begin());

        for(int idx = 0; idx != numpoints-2; ++idx) {
          real distance(boil::yotta);
          int newidx(0);
          for(int jdx = 0; jdx != pointset.size(); ++jdx) {
            XYZ distvector;
            distvector.x = pointset[jdx].x - orderedset[idx].x;
            distvector.y = pointset[jdx].y - orderedset[idx].y;
            distvector.z = pointset[jdx].z - orderedset[idx].z;

            real newdistance = distvector.x*distvector.x
                             + distvector.y*distvector.y
                             + distvector.z*distvector.z;

            newdistance = sqrt(newdistance);
            if(newdistance<distance) {
              distance = newdistance;
              newidx = jdx;
            }
          }
          orderedset.push_back(pointset[newidx]);
          pointset.erase(pointset.begin()+newidx);
        }
        /* last point is deterministic */
        orderedset.push_back(pointset[0]);

    
        /* calculate area */ 
        real polyarea = calc_area(orderedset,normvector); 

        adensgeom[i][j][k] = polyarea/dx/dy/dz;

#if 0
      boil::oout<<i<<" "<<j<<" "<<k<<" | "<<c;
      for(auto a: orderedset)
        boil::oout<<" | "<< a.x<<" "<<a.y<<" "<<a.z;
        boil::oout<<" | "<< adensgeom[i][j][k]<<" "<<adens[i][j][k];
      boil::oout<<boil::endl;
      boil::oout<<boil::endl;
#endif
      }
        
    }


  }

  real sum(0.0), sum2(0.0);
  int count(0), count2(0);
  for_vijk(adens,i,j,k) {
    real sumplus =  adens[i][j][k]*adens.dV(i,j,k);
    real sum2plus =  adensgeom[i][j][k]*adensgeom.dV(i,j,k);
    if(sumplus>boil::atto) count++;
    if(sum2plus>boil::atto) count2++;
    sum += sumplus;
    sum2 += sum2plus;
  }
  boil::oout<<"VOF::finescalar_adens "<<count<<" "<<sum<<boil::endl;
  boil::oout<<"VOF::finescalar_adensgeom "<<count2<<" "<<sum2<<boil::endl;

  boil::timer.stop("finescalar adens geom");

  return;
}
