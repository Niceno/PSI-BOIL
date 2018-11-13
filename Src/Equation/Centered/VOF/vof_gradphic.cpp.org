#include "vof.h"

/******************************************************************************/
void VOF::gradphic(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate grad(csa)/|grad(csa)| at cell center.
*         Resluts: nx, ny, nz
*******************************************************************************/

  /* cell centered base */
  for_ijk(i,j,k) {
    if(fabs(iflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
#if 0
      nx[i][j][k] = (sca[i+1][j][k]-sca[i-1][j][k])/(dxw(i)+dxe(i));
      ny[i][j][k] = (sca[i][j+1][k]-sca[i][j-1][k])/(dys(j)+dyn(j));
      nz[i][j][k] = (sca[i][j][k+1]-sca[i][j][k-1])/(dzb(k)+dzt(k));
#endif
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
      nx[i][j][k] = 0.25 * ( q100 - q000 + q110 - q010
                           + q101 - q001 + q111 - q011)/phi.dxc(i);
      ny[i][j][k] = 0.25 * ( q010 - q000 + q110 - q100
                           + q011 - q001 + q111 - q101)/phi.dyc(j);
      nz[i][j][k] = 0.25 * ( q001 - q000 + q101 - q100
                           + q011 - q010 + q111 - q110)/phi.dzc(k);
#endif
#if 0
      real vmm, vm, vp, vpp;
      // x-direction
      if (i==si()) { vmm=sca[i-1][j][k]; } else { vmm=sca[i-2][j][k]; }
      vm = sca[i-1][j][k];
      vp = sca[i+1][j][k];
      if (i==ei()) { vpp=sca[i+1][j][k]; } else { vpp=sca[i+2][j][k]; }
      nx[i][j][k] = (vmm-8.0*vm+8.0*vp-vmm)/(12.0*0.5*(dxw(i)+dxe(i)));

      // y-direction
      if (j==sj()) { vmm=sca[i][j-1][k]; } else { vmm=sca[i][j-2][k]; }
      vm = sca[i][j-1][k];
      vp = sca[i][j+1][k];
      if (j==ej()) { vpp=sca[i][j+1][k]; } else { vpp=sca[i][j+2][k]; }
      ny[i][j][k] = (vmm-8.0*vm+8.0*vp-vmm)/(12.0*0.5*(dys(j)+dyn(j)));

      // z-direction
      if (k==sk()) { vmm=sca[i][j][k-1]; } else { vmm=sca[i][j][k-2]; }
      vm = sca[i][j][k-1];
      vp = sca[i][j][k+1];
      if (k==ek()) { vpp=sca[i][j][k+1]; } else { vpp=sca[i][j][k+2]; }
      nz[i][j][k] = (vmm-8.0*vm+8.0*vp-vmm)/(12.0*0.5*(dzb(k)+dzt(k)));
#endif
    }
  }

  /* normal vector at adjacent cells next to wall, symmetric and IB */
  insert_bc_gradphic(sca); 

  /* normal vector on boundary plane */
  insert_bc_norm();

  /* normalize */
  for_avijk(sca,i,j,k) {
    if(fabs(iflag[i][j][k])>real(nlayer-4)){
      nx[i][j][k] = 0.0;
      ny[i][j][k] = 0.0;
      nz[i][j][k] = 0.0;
    } else {
      normalize(nx[i][j][k],ny[i][j][k],nz[i][j][k]);
    }
  }

  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

//  for (int k = 1; k < 32; k=k+1)
//  {
//    std::cout<<"gradphic "<<nx [16][32][k] <<"\n";
//  }
 
//  for_aijk(i,j,k)
//    std::cout<<"gradphic "<<nx [i][j][k] <<"\n";



//  boil::plot->plot(nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //boil::plot->plot(sca, "dist", time->current_step());
  //exit(0);

  return;
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_gradphic.cpp,v 1.7 2015/05/05 15:11:07 sato Exp $'/
+-----------------------------------------------------------------------------*/
