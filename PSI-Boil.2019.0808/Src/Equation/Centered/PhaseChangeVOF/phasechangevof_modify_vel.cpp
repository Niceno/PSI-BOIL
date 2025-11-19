#include "phasechangevof.h"
using namespace std;

/******************************************************************************/
void PhaseChangeVOF::modify_vel(Vector & uvw, 
                                const Vector & cnew, const Vector & cold) {
/***************************************************************************//**
*  \brief Modify velocity where interface crossed
    This function must be called after update() because M, nx, ny and nz
    must be computed before.
*******************************************************************************/
  //boil::oout<<"modify_vel:begin\n";
  Comp m = Comp::u();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=cold[m][i][j][k];
    real cf_new=cnew[m][i][j][k];
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        //std::cout<<"vel:cell1 "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
#if 0
        if (boil::cart.iam()==0) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
        }
        if (boil::cart.iam()==8) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
        }
#endif
        //std::cout<<"vel:cf   "<<boil::cart.iam()<<" "<<cf_new<<" "<<cf_old<<"\n";
        //std::cout<<"vel:nx   "<<boil::cart.iam()<<" "<<nz[i][j][k-1]<<" "<<nz[i][j][k]<<"\n";
        real nxf = 0.5*(nx[i-1][j][k]+nx[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        if (boil::cart.iam()==0) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
          std::cout<<"cold= "<<cold[i][j][k]<<" "<<cold[i-1][j][k]<<"\n";
          std::cout<<"M= "<<M[i][j][k]<<" "<<M[i-1][j][k]<<"\n";
        }
        if (boil::cart.iam()==8) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
          std::cout<<"cold= "<<cold[i][j][k]<<" "<<cold[i-1][j][k]<<"\n";
          std::cout<<"M= "<<M[i][j][k]<<" "<<M[i-1][j][k]<<"\n";
        }
#endif
#if 0
        real Mdot;
        if ((0.5-cold[i-1][j][k])/(cold[i][j][k]-cold[i-1][j][k])<0.5) {
          Mdot=M[i-1][j][k];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i-1][j][k]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nxf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nxf*(1.0/rhof_old-1.0/rhof_new);
        //std::cout<<"vel:cell2 "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
  }

  //boil::oout<<"modify_vel:v\n";
  m = Comp::v();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=cold[m][i][j][k];
    real cf_new=cnew[m][i][j][k];
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        real nyf = 0.5*(ny[i][j-1][k]+ny[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        real Mdot;
        if ((0.5-cold[i][j-1][k])/(cold[i][j][k]-cold[i][j-1][k])<0.5) {
          Mdot=M[i][j-1][k];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i][j-1][k]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nyf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nyf*(1.0/rhof_old-1.0/rhof_new);
      }
    }
  }

  //boil::oout<<"modify_vel:w\n";
  m = Comp::w();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=cold[m][i][j][k];
    real cf_new=cnew[m][i][j][k];
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        //std::cout<<"vel:cell "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
        //std::cout<<"vel:cf   "<<boil::cart.iam()<<" "<<cf_new<<" "<<cf_old<<"\n";
        //std::cout<<"vel:nx   "<<boil::cart.iam()<<" "<<nz[i][j][k-1]<<" "<<nz[i][j][k]<<"\n";
        real nzf = 0.5*(nz[i][j][k-1]+nz[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        real Mdot;
        if ((0.5-cold[i][j][k-1])/(cold[i][j][k]-cold[i][j][k-1])<0.5) {
          Mdot=M[i][j][k-1];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i][j][k-1]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nzf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nzf*(1.0/rhof_old-1.0/rhof_new);
      }
    }
  }

  //boil::oout<<"modify_vel:exchange\n";
  //uvw->exchange();
  uvw.exchange();
  //boil::oout<<"modify_vel:end\n";

}

/******************************************************************************/
void PhaseChangeVOF::modify_vel(Vector & uvw, 
                                const Scalar & cnew, const Scalar & cold) {
/***************************************************************************//**
*  \brief Modify velocity where interface crossed
    This function must be called after update() because M, nx, ny and nz
    must be computed before.
*******************************************************************************/
  //boil::oout<<"modify_vel:begin\n";
  Comp m = Comp::u();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=0.5*(cold[i-1][j][k]+cold[i][j][k]);
    real cf_new=0.5*(cnew[i-1][j][k]+cnew[i][j][k]);
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        //std::cout<<"vel:cell1 "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
#if 0
        if (boil::cart.iam()==0) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
        }
        if (boil::cart.iam()==8) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
        }
#endif
        //std::cout<<"vel:cf   "<<boil::cart.iam()<<" "<<cf_new<<" "<<cf_old<<"\n";
        //std::cout<<"vel:nx   "<<boil::cart.iam()<<" "<<nz[i][j][k-1]<<" "<<nz[i][j][k]<<"\n";
        real nxf = 0.5*(nx[i-1][j][k]+nx[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        if (boil::cart.iam()==0) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
          std::cout<<"cold= "<<cold[i][j][k]<<" "<<cold[i-1][j][k]<<"\n";
          std::cout<<"M= "<<M[i][j][k]<<" "<<M[i-1][j][k]<<"\n";
        }
        if (boil::cart.iam()==8) {
          std::cout<<"nx[i  ]= "<<boil::cart.iam()<<" "<<i<<" "<<nx[i][j][k]<<"\n";
          std::cout<<"nx[i-1]= "<<boil::cart.iam()<<" "<<i-1<<" "<<nx[i-1][j][k]<<"\n";
          std::cout<<"cf_new= "<<boil::cart.iam()<<" "<<cf_new<<"\n";
          std::cout<<"cold= "<<cold[i][j][k]<<" "<<cold[i-1][j][k]<<"\n";
          std::cout<<"M= "<<M[i][j][k]<<" "<<M[i-1][j][k]<<"\n";
        }
#endif
#if 0
        real Mdot;
        if ((0.5-cold[i-1][j][k])/(cold[i][j][k]-cold[i-1][j][k])<0.5) {
          Mdot=M[i-1][j][k];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i-1][j][k]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nxf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nxf*(1.0/rhof_old-1.0/rhof_new);
        //std::cout<<"vel:cell2 "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
  }

  //boil::oout<<"modify_vel:v\n";
  m = Comp::v();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=0.5*(cold[i][j-1][k]+cold[i][j][k]);
    real cf_new=0.5*(cnew[i][j-1][k]+cnew[i][j][k]);
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        real nyf = 0.5*(ny[i][j-1][k]+ny[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        real Mdot;
        if ((0.5-cold[i][j-1][k])/(cold[i][j][k]-cold[i][j-1][k])<0.5) {
          Mdot=M[i][j-1][k];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i][j-1][k]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nyf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nyf*(1.0/rhof_old-1.0/rhof_new);
      }
    }
  }

  //boil::oout<<"modify_vel:w\n";
  m = Comp::w();
  //for_vmijk((*uvw), m,i,j,k) {
  for_vmijk(uvw, m,i,j,k) {
    real cf_old=0.5*(cold[i][j][k-1]+cold[i][j][k]);
    real cf_new=0.5*(cnew[i][j][k-1]+cnew[i][j][k]);
    if ((cf_old-clrsurf)*(cf_new-clrsurf)<0.0) {
      if (dom->ibody().on(m,i,j,k)) {
        //std::cout<<"vel:cell "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
        //std::cout<<"vel:cf   "<<boil::cart.iam()<<" "<<cf_new<<" "<<cf_old<<"\n";
        //std::cout<<"vel:nx   "<<boil::cart.iam()<<" "<<nz[i][j][k-1]<<" "<<nz[i][j][k]<<"\n";
        real nzf = 0.5*(nz[i][j][k-1]+nz[i][j][k]);
        real rhof_new, rhof_old;
        if(cf_new<0.5){
          rhof_new=rhov;
          rhof_old=rhol;
        } else {
          rhof_new=rhol;
          rhof_old=rhov;
        }
#if 0
        real Mdot;
        if ((0.5-cold[i][j][k-1])/(cold[i][j][k]-cold[i][j][k-1])<0.5) {
          Mdot=M[i][j][k-1];
        } else {
          Mdot=M[i][j][k];
        }
#else
        real Mdot = 0.5*(M[i][j][k-1]+M[i][j][k]);
#endif
        //(*uvw)[m][i][j][k] += Mdot*nzf*(1.0/rhof_old-1.0/rhof_new);
        uvw[m][i][j][k] += Mdot*nzf*(1.0/rhof_old-1.0/rhof_new);
      }
    }
  }

  //boil::oout<<"modify_vel:exchange\n";
  //uvw->exchange();
  uvw.exchange();
  //boil::oout<<"modify_vel:end\n";

}

