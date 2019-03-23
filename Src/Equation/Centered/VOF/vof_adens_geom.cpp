#include "vof.h"

/* calculate the cross product of the arguments */
VOF::XYZ VOF::CrossProduct(const XYZ p1, const XYZ p2) {
  XYZ cp;
  cp.x = p1.y * p2.z - p1.z * p2.y;
  cp.y = p1.z * p2.x - p1.x * p2.z;
  cp.z = p1.x * p2.y - p1.y * p2.x;

  return(cp);
}

/* calculate the dot product of the arguments */
real VOF::DotProduct(const XYZ p1, const XYZ p2) {
  real dp(0.0);
  dp += p1.x * p2.x;
  dp += p1.y * p2.y;
  dp += p1.z * p2.z;

  return(dp);
}

/* add two xyz vectors */
VOF::XYZ VOF::PlusXYZ(const XYZ p1, const XYZ p2) {
  XYZ pv;
  pv.x = p1.x + p2.x;
  pv.y = p1.y + p2.y;
  pv.z = p1.z + p2.z;

  return(pv);
}

real VOF::calc_area(const std::vector<XYZ> &vect, const XYZ norm) { 
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
void VOF::cal_adens_geom(Scalar & eval) {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adensgeom
*******************************************************************************/

#if 0
  /* debugging */
  stmp=0.0;
  stmp2=0.0;
#endif
 

  /* cell centered */
  for_vijk(eval,i,j,k) {
    bool interface = eval[i][j][k]>boil::pico;
    bool vicinity  = eval[i+1][j][k]>boil::pico;
    vicinity |= eval[i-1][j][k]>boil::pico;
    vicinity |= eval[i][j+1][k]>boil::pico;
    vicinity |= eval[i][j-1][k]>boil::pico;
    vicinity |= eval[i][j][k+1]>boil::pico;
    vicinity |= eval[i][j][k-1]>boil::pico;
 
    if(interface||(vicinity&&phi[i][j][k]-1.<-boil::micro)) {
      /* n points to the gas, in normalized space */
      real n1 = -nx[i][j][k];
      real n2 = -ny[i][j][k];
      real n3 = -nz[i][j][k];
      /* m points to the gas, in real space */
      real m1 = -mx[i][j][k];
      real m2 = -my[i][j][k];
      real m3 = -mz[i][j][k];

      /* standardized space alpha */
      real alpha = nalpha[i][j][k];

      /* vn is positive and standardized */
      real vn1 = fabs(n1);
      real vn2 = fabs(n2);
      real vn3 = fabs(n3);

      /* vm is positive */
      real vm1 = fabs(m1);
      real vm2 = fabs(m2);
      real vm3 = fabs(m3);

#if 1
      real dx, dy, dz;

      real vnmin, vmmin;
      real vnmid, vmmid;
      real vnmax, vmmax;
      if       (vn1<=vn2&&vn1<=vn3) {
        vnmin = vn1; 
        vmmin = vm1;
        dz = eval.dxc(i);
        if(vn2<=vn3) {
          vnmax = vn3;
          vmmax = vm3;
          dx = eval.dzc(k);
          vnmid = vn2;
          vmmid = vm2;
          dy = eval.dyc(j);
        } else {
          vnmax = vn2;
          vmmax = vm2;
          dx = eval.dyc(j);
          vnmid = vn3;
          vmmid = vm3;
          dy = eval.dzc(k);
        }
      } else if(vn2<=vn1&&vn2<=vn3) {
        vnmin = vn2; 
        vmmin = vm2;
        dz = eval.dyc(j);
        if(vn1<=vn3) {
          vnmax = vn3;
          vmmax = vm3;
          dx = eval.dzc(k);
          vnmid = vn1;
          vmmid = vm1;
          dy = eval.dxc(i);
        } else {
          vnmax = vn1;
          vmmax = vm1;
          dx = eval.dxc(i);
          vnmid = vn3;
          vmmid = vm3;
          dy = eval.dzc(k);
        }
      } else {
        vnmin = vn3; 
        vmmin = vm3;
        dz = eval.dzc(k);
        if(vn1<=vn2) {
          vnmax = vn2;
          vmmax = vm2;
          dx = eval.dyc(j);
          vnmid = vn1;
          vmmid = vm1;
          dy = eval.dxc(i);
        } else {
          vnmax = vn1;
          vmmax = vm1;
          dx = eval.dxc(i);
          vnmid = vn2;
          vmmid = vm2;
          dy = eval.dyc(j);
        }
      }

      vn3 = vnmin;
      vn2 = vnmid;
      vn1 = vnmax;
      vm3 = vmmin;
      vm2 = vmmid;
      vm1 = vmmax;
#endif

      /* beyond this point, x and y and z don't have their usual meanings:
       * rather x corresponds to the largest normal vector components and so on */

      XYZ normvector;
      normvector.x = vm1;
      normvector.y = vm2;
      normvector.z = vm3;

      /* for each edge of the cell, try to find an intersection */
      std::list<XYZ> pointset;

      bool xnorm_nonzero(vn1>boil::pico);
      bool ynorm_nonzero(vn2>boil::pico);
      bool znorm_nonzero(vn3>boil::pico);
      
      real xval, yval, zval;

      /* y = 0, z = 0, x varies */
      yval = 0.0;
      zval = 0.0;
 
      if(xnorm_nonzero) { 
        xval = (alpha-vn2*yval-vn3*zval)/vn1;
        if(xval>=0.0&&xval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval; 
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        } 
      }

      /* y = 0, z = dz, x varies */
      yval = 0.0;
      zval = 1.0;

      if(xnorm_nonzero) {
        xval = (alpha-vn2*yval-vn3*zval)/vn1;
        if(xval>=0.0&&xval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* y = dy, z = 0, x varies */
      yval = 1.0;
      zval = 0.0;

      if(xnorm_nonzero) {
        xval = (alpha-vn2*yval-vn3*zval)/vn1;
        if(xval>=0.0&&xval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* y = dy, z = dz, x varies */
      yval = 1.0;
      zval = 1.0;

      if(xnorm_nonzero) {
        xval = (alpha-vn2*yval-vn3*zval)/vn1;
        if(xval>=0.0&&xval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, z = 0, y varies */
      xval = 0.0;
      zval = 0.0;

      if(ynorm_nonzero) {
        yval = (alpha-vn1*xval-vn3*zval)/vn2;
        if(yval>=0.0&&yval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, z = dz, y varies */
      xval = 0.0;
      zval = 1.0;

      if(ynorm_nonzero) {
        yval = (alpha-vn1*xval-vn3*zval)/vn2;
        if(yval>=0.0&&yval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, z = 0, y varies */
      xval = 1.0;
      zval = 0.0;

      if(ynorm_nonzero) {
        yval = (alpha-vn1*xval-vn3*zval)/vn2;
        if(yval>=0.0&&yval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, z = dz, y varies */
      xval = 1.0;
      zval = 1.0;

      if(ynorm_nonzero) {
        yval = (alpha-vn1*xval-vn3*zval)/vn2;
        if(yval>=0.0&&yval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, y = 0, z varies */
      xval = 0.0;
      yval = 0.0;

      if(znorm_nonzero) {
        zval = (alpha-vn1*xval-vn2*yval)/vn3;
        if(zval>=0.0&&zval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, y = 0, z varies */
      xval = 1.0;
      yval = 0.0;

      if(znorm_nonzero) {
        zval = (alpha-vn1*xval-vn2*yval)/vn3;
        if(zval>=0.0&&zval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = 0, y = dy, z varies */
      xval = 0.0;
      yval = 1.0;

      if(znorm_nonzero) {
        zval = (alpha-vn1*xval-vn2*yval)/vn3;
        if(zval>=0.0&&zval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

      /* x = dx, y = dy, z varies */
      xval = 1.0;
      yval = 1.0;

      if(znorm_nonzero) {
        zval = (alpha-vn1*xval-vn2*yval)/vn3;
        if(zval>=0.0&&zval<=1.0) {
          XYZ newpoint;
          newpoint.x = dx*xval;
          newpoint.y = dy*yval;
          newpoint.z = dz*zval;
          pointset.push_back(newpoint);
        }
      }

#if 0
      /* debugging */
      if(fabs(adens[i][j][k]-195111.)<1.0) {
        boil::oout<<int(pointset.size())<<" "<<vn1<<" "<<vn2<<" "<<vn3<<" | "<<n1<<" "<<n2<<" "<<n3;
        for(std::list<XYZ>::iterator p = pointset.begin();
            p != pointset.end(); ++p) {
          boil::oout<<" | "<<(*p).x<<" "<<(*p).y<<" "<<(*p).z;
        }
        boil::oout<<boil::endl;
        boil::oout<<boil::endl;
      }
#endif

#if 0
      /* debugging */
      stmp[i][j][k] = pointset.size();
      stmp2[i][j][k] = alpha;  
#endif

      if(pointset.size()<3) {
        adensgeom[i][j][k] = 0.0;
        //boil::oout<<"FS::caladensgeom: Warning, inconsistent geometry"<<boil::endl;
      } else {
        /* order points */
        /* for a convex polygon, neighboring vertices are the closest ones */
        std::vector<XYZ> orderedset;
        int numpoints = pointset.size();

        /* first point is arbitrary */ 
        std::list<XYZ>::iterator b = pointset.begin(); 
        orderedset.push_back(*b);
        pointset.erase(b);

        for(int idx = 0; idx != numpoints-2; ++idx) {
          real distance(boil::yotta);
          std::list<XYZ>::iterator n;
          int newidx(0);
          for(std::list<XYZ>::iterator p = pointset.begin();
              p != pointset.end(); ++p) {
            XYZ distvector;
            distvector.x = (*p).x - orderedset[idx].x;
            distvector.y = (*p).y - orderedset[idx].y;
            distvector.z = (*p).z - orderedset[idx].z;

            real newdistance = distvector.x*distvector.x
                             + distvector.y*distvector.y
                             + distvector.z*distvector.z;

            if(newdistance<distance) {
              distance = newdistance;
              n = p;
            }
          }
          orderedset.push_back(*n);
          pointset.erase(n);
        }
        /* last point is deterministic */
        orderedset.push_back(*(pointset.begin()));

    
        /* calculate area */ 
        real polyarea = calc_area(orderedset,normvector); 

        adensgeom[i][j][k] = polyarea/dx/dy/dz;

#if 0
      if(  ((i==40) && (j==40) && (k==53))
         ||((i==44) && (j==40) && (k==40))) {
        boil::oout<<i<<" "<<j<<" "<<k<<" | "<<phi[i][j][k]<<" | "<<vn1<<" "<<vn2<<" "<<vn3<<" | "<<nalpha[i][j][k];
        for(auto a: orderedset)
          boil::oout<<" | "<< a.x<<" "<<a.y<<" "<<a.z;
        boil::oout<<" | "<< adensgeom[i][j][k]<<" "<<adens[i][j][k];
        boil::oout<<boil::endl;
        boil::oout<<boil::endl;
      }
#endif
      }
        
    } /* eval nonzero */
  } /* for vijk */
  adensgeom.exchange();

  real sum(0.0), sum2(0.0);
  int count(0), count2(0);
  for_vijk(eval,i,j,k) {
    real sumplus =  eval[i][j][k]*eval.dV(i,j,k);
    real sum2plus =  adensgeom[i][j][k]*adensgeom.dV(i,j,k);
    if(sumplus>boil::atto) count++;
    if(sum2plus>boil::atto) count2++;
    sum += sumplus;
    sum2 += sum2plus;
  }
  boil::oout<<"VOF::finescalar_adens "<<count<<" "<<sum<<boil::endl;
  boil::oout<<"VOF::finescalar_adensgeom "<<count2<<" "<<sum2<<boil::endl;

#if 1
  boil::plot->plot(phi,eval,adensgeom,"clr-adens-adensgeom",time->current_step());
  //exit(0);
#endif


  return;
}
