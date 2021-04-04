#include "vof.h"
#include <numeric>

//#define USE_VICINITY
//#define USE_TOLERANCE

struct XYZ {
  real x,y,z;
};

/* calculate the cross product of the arguments */
XYZ CrossProduct(const XYZ p1, const XYZ p2) {
  XYZ cp;
  cp.x = p1.y * p2.z - p1.z * p2.y;
  cp.y = p1.z * p2.x - p1.x * p2.z;
  cp.z = p1.x * p2.y - p1.y * p2.x;

  return(cp);
}

/* calculate the dot product of the arguments */
real DotProduct(const XYZ p1, const XYZ p2) {
  real dp(0.0);
  dp += p1.x * p2.x;
  dp += p1.y * p2.y;
  dp += p1.z * p2.z;

  return(dp);
}

/* add two xyz vectors */
XYZ PlusXYZ(const XYZ p1, const XYZ p2) {
  XYZ pv;
  pv.x = p1.x + p2.x;
  pv.y = p1.y + p2.y;
  pv.z = p1.z + p2.z;

  return(pv);
}

/* subtract two xyz vectors */
XYZ MinusXYZ(const XYZ p1, const XYZ p2) {
  XYZ pv;
  pv.x = p1.x - p2.x;
  pv.y = p1.y - p2.y;
  pv.z = p1.z - p2.z;

  return(pv);
}

XYZ Normalize(const XYZ p) {
  XYZ n;
  real norm = p.x*p.x + p.y*p.y + p.z*p.z;
  norm = sqrt(norm);
  n.x = p.x/norm;
  n.y = p.y/norm;
  n.z = p.z/norm;

  return n;
}

real calc_area(const std::vector<XYZ> &vect, const XYZ norm) { 
  XYZ cross;
  cross.x = 0.0;
  cross.y = 0.0;
  cross.z = 0.0;
  for(int i = 0; i != vect.size(); ++i) {
    cross = PlusXYZ(cross,CrossProduct(vect[i],vect[(i+1) % vect.size()]));
  }
  return 0.5 * fabs(DotProduct(cross,norm));
}

/* https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-
 *       in-the-same-way-with-criteria-that-uses-only-one-of */
std::vector<std::size_t> sort_permutation(const std::vector<real>& vec) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
      [&](std::size_t i, std::size_t j){ return vec[i]<vec[j]; });
  return p;
}

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}

/******************************************************************************/
void VOF::cal_adens_geom(Scalar & adensgeom, const Scalar & sca,
                         const bool use_vicinity) {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adensgeom
*
*  Used only for comparison purposes !!!!!!!!!
*******************************************************************************/

#if 0
  /* debugging */
  stmp=0.0;
  stmp2=0.0;
#endif

  adensgeom = 0.0;

  /* cell centered */
  for_vijk(adensgeom,i,j,k) {
    bool condition(true);
    real valcc = sca[i  ][j  ][k  ];
    real valmc = sca[i-1][j  ][k  ];
    real valpc = sca[i+1][j  ][k  ];
    real valcm = sca[i  ][j-1][k  ];
    real valcp = sca[i  ][j+1][k  ];
    real valmm = sca[i  ][j  ][k-1];
    real valpm = sca[i  ][j  ][k+1];

    real val1 = sca[i-1][j-1][k  ];
    real val2 = sca[i-1][j+1][k  ];
    real val3 = sca[i-1][j  ][k-1];
    real val4 = sca[i-1][j  ][k+1];
    real val5 = sca[i+1][j-1][k  ];
    real val6 = sca[i+1][j+1][k  ];
    real val7 = sca[i+1][j  ][k-1];
    real val8 = sca[i+1][j  ][k+1];

    real val9 = sca[i  ][j-1][k-1];
    real valA = sca[i  ][j-1][k+1];
    real valB = sca[i  ][j+1][k-1];
    real valC = sca[i  ][j+1][k+1];

    real valD = sca[i+1][j-1][k-1];
    real valE = sca[i+1][j-1][k+1];
    real valF = sca[i+1][j+1][k-1];
    real valG = sca[i+1][j+1][k+1];

    real valH = sca[i-1][j-1][k-1];
    real valI = sca[i-1][j-1][k+1];
    real valJ = sca[i-1][j+1][k-1];
    real valK = sca[i-1][j+1][k+1];

    bool vicinity = ((valcc-phisurf)*(valmc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcp-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valmm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpm-phisurf)<=0.0)
#if 1
                  ||((valcc-phisurf)*(val1-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val2-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val3-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val4-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val5-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val6-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val7-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val8-phisurf)<=0.0)
                  ||((valcc-phisurf)*(val9-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valA-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valB-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valC-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valD-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valE-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valF-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valG-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valH-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valI-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valJ-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valK-phisurf)<=0.0)
#endif
                  ;

    if(use_vicinity)
      condition &= vicinity;
#ifdef USE_TOLERANCE
    condition &= (valcc-1.<-boil::micro)&&valcc>boil::micro;
#endif
    if(condition) {
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

      if(!boil::realistic(alpha))
        continue;

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
        dz = phi.dxc(i);
        if(vn2<=vn3) {
          vnmax = vn3;
          vmmax = vm3;
          dx = phi.dzc(k);
          vnmid = vn2;
          vmmid = vm2;
          dy = phi.dyc(j);
        } else {
          vnmax = vn2;
          vmmax = vm2;
          dx = phi.dyc(j);
          vnmid = vn3;
          vmmid = vm3;
          dy = phi.dzc(k);
        }
      } else if(vn2<=vn1&&vn2<=vn3) {
        vnmin = vn2; 
        vmmin = vm2;
        dz = phi.dyc(j);
        if(vn1<=vn3) {
          vnmax = vn3;
          vmmax = vm3;
          dx = phi.dzc(k);
          vnmid = vn1;
          vmmid = vm1;
          dy = phi.dxc(i);
        } else {
          vnmax = vn1;
          vmmax = vm1;
          dx = phi.dxc(i);
          vnmid = vn3;
          vmmid = vm3;
          dy = phi.dzc(k);
        }
      } else {
        vnmin = vn3; 
        vmmin = vm3;
        dz = phi.dzc(k);
        if(vn1<=vn2) {
          vnmax = vn2;
          vmmax = vm2;
          dx = phi.dyc(j);
          vnmid = vn1;
          vmmid = vm1;
          dy = phi.dxc(i);
        } else {
          vnmax = vn1;
          vmmax = vm1;
          dx = phi.dxc(i);
          vnmid = vn2;
          vmmid = vm2;
          dy = phi.dyc(j);
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
      std::vector<XYZ> pointset;

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
        for(auto p : pointset) {
          boil::oout<<" | "<<p.x<<" "<<p.y<<" "<<p.z;
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
        //boil::oout<<"FS::caladensgeom: Warning, inconsistent geometry "<<alpha<<boil::endl;
      } else {
        /* order points */
        std::vector<XYZ> orderedset;
        int numpoints = pointset.size();

        /* find centroid */
        XYZ centroid;
        centroid.x = centroid.y = centroid.z = 0.;
        for(auto & p : pointset) {
          centroid.x += p.x/numpoints;
          centroid.y += p.y/numpoints;
          centroid.z += p.z/numpoints;
        }

        /* first point is arbitrary */ 
        auto b = pointset.begin(); 
        orderedset.push_back(*b);
        pointset.erase(b);

        /* vector to first point */
        XYZ vectfirst = MinusXYZ(orderedset[0],centroid);

        /* angles of remaining points */
        std::vector<real> angles;

        for(auto & p : pointset) {
          XYZ vectp = MinusXYZ(p,centroid);
          real angle = acos(std::min(1.,std::max(-1.,
                                     DotProduct(Normalize(vectfirst),
                                                Normalize(vectp))
                           )));
          XYZ cross = CrossProduct(vectfirst,vectp);
          if(DotProduct(normvector,cross)<0.) {
            angle = 2*boil::pi-angle;
          }
          angles.push_back(angle);
        }

        /* sort according to angles */
        auto per = sort_permutation(angles);
        pointset = apply_permutation(pointset,per);

        /* add to orderedset */
        orderedset.insert(orderedset.end(),pointset.begin(),pointset.end());

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
        
    } /* condition true */
  } /* for vijk */
  adensgeom.exchange();

  return;
}

/***********************************************************************/
real local_gradclr(const Scalar & sca, const int i, const int j, const int k) {
/***********************************************************************/

#if 1
  real q000, q001, q010, q011, q100, q101, q110, q111;
  q000 = 0.125 * (sca[i-1][j-1][k-1] + sca[i][j-1][k-1]
                + sca[i-1][j  ][k-1] + sca[i][j  ][k-1]
                + sca[i-1][j-1][k  ] + sca[i][j-1][k  ]
                + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]);
  q001 = 0.125 * (sca[i-1][j-1][k  ] + sca[i][j-1][k  ]
                + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                + sca[i-1][j-1][k+1] + sca[i][j-1][k+1]
                + sca[i-1][j  ][k+1] + sca[i][j  ][k+1]);
  q100 = 0.125 * (sca[i  ][j-1][k-1] + sca[i+1][j-1][k-1]
                + sca[i  ][j  ][k-1] + sca[i+1][j  ][k-1]
                + sca[i  ][j-1][k  ] + sca[i+1][j-1][k  ]
                + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]);
  q101 = 0.125 * (sca[i  ][j-1][k  ] + sca[i+1][j-1][k  ]
                + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                + sca[i  ][j-1][k+1] + sca[i+1][j-1][k+1]
                + sca[i  ][j  ][k+1] + sca[i+1][j  ][k+1]);
  q010 = 0.125 * (sca[i-1][j  ][k-1] + sca[i][j  ][k-1]
                + sca[i-1][j+1][k-1] + sca[i][j+1][k-1]
                + sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                + sca[i-1][j+1][k  ] + sca[i][j+1][k  ]);
  q011 = 0.125 * (sca[i-1][j  ][k  ] + sca[i][j  ][k  ]
                + sca[i-1][j+1][k  ] + sca[i][j+1][k  ]
                + sca[i-1][j  ][k+1] + sca[i][j  ][k+1]
                + sca[i-1][j+1][k+1] + sca[i][j+1][k+1]);
  q110 = 0.125 * (sca[i  ][j  ][k-1] + sca[i+1][j  ][k-1]
                + sca[i  ][j+1][k-1] + sca[i+1][j+1][k-1]
                + sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                + sca[i  ][j+1][k  ] + sca[i+1][j+1][k  ]);
  q111 = 0.125 * (sca[i  ][j  ][k  ] + sca[i+1][j  ][k  ]
                + sca[i  ][j+1][k  ] + sca[i+1][j+1][k  ]
                + sca[i  ][j  ][k+1] + sca[i+1][j  ][k+1]
                + sca[i  ][j+1][k+1] + sca[i+1][j+1][k+1]);

  real gradx = 0.25 * ( q100 - q000 + q110 - q010
                       + q101 - q001 + q111 - q011)/sca.dxc(i);
  real grady = 0.25 * ( q010 - q000 + q110 - q100
                       + q011 - q001 + q111 - q101)/sca.dyc(j);
  real gradz = 0.25 * ( q001 - q000 + q101 - q100
                       + q011 - q010 + q111 - q110)/sca.dzc(k);
#else
  real gradx = (sca[i+1][j][k]-sca[i-1][j][k])/(sca.dxw(i)+sca.dxe(i));
  real grady = (sca[i][j+1][k]-sca[i][j-1][k])/(sca.dys(j)+sca.dyn(j));
  real gradz = (sca[i][j][k+1]-sca[i][j][k-1])/(sca.dzb(k)+sca.dzt(k));
#endif

  return sqrt(gradx*gradx+grady*grady+gradz*gradz);
}

/***********************************************************************/
void VOF::cal_adens_gradclr(Scalar & adensgeom, const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adensgeom
*
*  Used only for comparison purposes !!!!!!!!!
*******************************************************************************/

  adensgeom = 0.0;

  /* cell centered */
  for_vijk(adensgeom,i,j,k) {
    bool condition(true);
#ifdef USE_VICINITY
    real valcc = sca[i  ][j  ][k  ];
    real valmc = sca[i-1][j  ][k  ];
    real valpc = sca[i+1][j  ][k  ];
    real valcm = sca[i  ][j-1][k  ];
    real valcp = sca[i  ][j+1][k  ];
    real valmm = sca[i  ][j  ][k-1];
    real valpm = sca[i  ][j  ][k+1];
    bool vicinity = ((valcc-phisurf)*(valmc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcp-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valmm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpm-phisurf)<=0.0);

    condition &= vicinity;
#endif
#ifdef USE_TOLERANCE
    condition &= (valcc-1.<-boil::micro)&&valcc>boil::micro;
#endif
    if(condition) {
      adensgeom[i][j][k] = local_gradclr(sca,i,j,k);
    }
  }

  return;
}

/***********************************************************************/
void VOF::cal_adens_gradclr_2phi(Scalar & adensgeom, const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adensgeom
*
*  Used only for comparison purposes !!!!!!!!!
*******************************************************************************/

  adensgeom = 0.0;

  /* cell centered */
  for_vijk(adensgeom,i,j,k) {
    bool condition(true);
#ifdef USE_VICINITY
    real valcc = sca[i  ][j  ][k  ];
    real valmc = sca[i-1][j  ][k  ];
    real valpc = sca[i+1][j  ][k  ];
    real valcm = sca[i  ][j-1][k  ];
    real valcp = sca[i  ][j+1][k  ];
    real valmm = sca[i  ][j  ][k-1];
    real valpm = sca[i  ][j  ][k+1];
    bool vicinity = ((valcc-phisurf)*(valmc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcp-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valmm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpm-phisurf)<=0.0);

    condition &= vicinity;
#endif
#ifdef USE_TOLERANCE
    condition &= (valcc-1.<-boil::micro)&&valcc>boil::micro;
#endif
    if(condition) {
      adensgeom[i][j][k] = 2.*std::min(1.0,std::max(0.0,sca[i][j][k]))
                  *local_gradclr(sca,i,j,k);
    }
  }

  return;
}

/***********************************************************************/
void VOF::cal_adens_gradclr_6phi(Scalar & adensgeom, const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate |grad(clr)| at cell center.
*         Results: gradclr = adensgeom
*
*  Used only for comparison purposes !!!!!!!!!
*******************************************************************************/

  adensgeom = 0.0;

  /* cell centered */
  for_vijk(adensgeom,i,j,k) {
    bool condition(true);
#ifdef USE_VICINITY
    real valcc = sca[i  ][j  ][k  ];
    real valmc = sca[i-1][j  ][k  ];
    real valpc = sca[i+1][j  ][k  ];
    real valcm = sca[i  ][j-1][k  ];
    real valcp = sca[i  ][j+1][k  ];
    real valmm = sca[i  ][j  ][k-1];
    real valpm = sca[i  ][j  ][k+1];
    bool vicinity = ((valcc-phisurf)*(valmc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpc-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valcp-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valmm-phisurf)<=0.0)
                  ||((valcc-phisurf)*(valpm-phisurf)<=0.0);

    condition &= vicinity;
#endif
#ifdef USE_TOLERANCE
    condition &= (valcc-1.<-boil::micro)&&valcc>boil::micro;
#endif
    if(condition) {
      real sc = std::min(1.0,std::max(0.0,sca[i][j][k]));
      adensgeom[i][j][k] = 6.*sc*(1.-sc)*local_gradclr(sca,i,j,k);
    }
  }

  return;
}
