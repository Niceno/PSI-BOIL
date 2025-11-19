#include "body.h"

/******************************************************************************/
real area(const real x0, const real y0, 
          const real x1, const real y1,
          const real x2, const real y2) {

  return 0.5 * (  ( y0 + y1 ) * ( x0 - x1 ) +
                  ( y1 + y2 ) * ( x1 - x2 ) +
                  ( y2 + y0 ) * ( x2 - x0 ) ) ;
}

/******************************************************************************/
void arrange_poly(const int nn, 
                  real x[], real y[], real z[],
                  real nx, real ny, real nz) {

   /* normalize */
   const real nmag=sqrt(nx*nx+ny*ny+nz*nz);
   nx=nx/nmag;
   ny=ny/nmag;
   nz=nz/nmag;

   /*-----------------------------------------+ 
   |  take absolute values of surface normal  |
   +-----------------------------------------*/
   const real anx = fabs(nx);
   const real any = fabs(ny);
   const real anz = fabs(nz);
   const real maxan = boil::maxr( anx, any, anz );

   /*--------------------+
   |  center of gravity  |
   +--------------------*/
   real g[3] = {0.0, 0.0, 0.0}; /* gx, gy, gz */
   for(int n=0; n<nn; n++) {
     g[0] += x[n];
     g[1] += y[n];
     g[2] += z[n];
   }
   for(int i=0; i<3; i++) g[i] /= (real)nn;
     
   std::vector<AngleNode> nodes;

   /*-----------------------+
   |  sort nodes by angles  |
   +-----------------------*/
     
   // project to plane, which pass center of gravity
   real d = -(nx*g[0] + ny*g[1] + nz*g[2]);
   real len0 = nx*x[0] + ny*y[0] + nz*z[0] + d;
   real xpro0= x[0] - len0*nx;
   real ypro0= y[0] - len0*ny;
   real zpro0= z[0] - len0*nz;

#if 0
   // check projected point is on the plane
   std::cout<<"on plane? "<< nx*xpro0+ny*ypro0+nz*zpro0+d<<"\n";
#endif

   // element vector e1
   real abse1=sqrt((xpro0-g[0])*(xpro0-g[0])
                  +(ypro0-g[1])*(ypro0-g[1])
                  +(zpro0-g[2])*(zpro0-g[2]));
   assert(abse1!=0.0);
   real e1x=(xpro0-g[0])/abse1;
   real e1y=(ypro0-g[1])/abse1;
   real e1z=(zpro0-g[2])/abse1;
   real e2x, e2y, e2z;
   // element vector e2
   if(e1y*nz-e1z*ny!=0.0){
     e2x=1.0;
     e2y=(nz*e1x*e2x-nx*e1z*e2x)/(ny*e1z-e1y*nz);
     e2z=(ny*e1x*e2x-nx*e1y*e2x)/(e1y*nz-ny*e1z);
   } else if(nz*e1x-nx*e1z!=0){
     e2y=1.0;
     e2x=(nz*e1y*e2y-ny*e1z*e2y)/(nx*e1z-nz*e1x);
     e2z=(nx*e1y*e2y-ny*e1x*e2y)/(nz*e1x-nx*e1z);
   } else if(ny*e1x-nx*e1y!=0){
     e2z=1.0;
     e2x=(ny*e1z*e2z-nz*e1y*e2z)/(nx*e1y-ny*e1x);
     e2y=(nx*e1z*e2z-nz*e1x*e2z)/(ny*e1x-nx*e1y);
   } else {
     std::cout<<"arrange_poly: Error!!!\n";
     std::cout<<"nodes lie on a line.\n";
     std::cout<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<"\n";
     std::cout<<"d= "<<d<<" len0= "<<len0<<"\n";
     std::cout<<"nn= "<<nn<<"\n";
     std::cout<<"p0 "<<x[0]<<" "<<y[0]<<" "<<z[0]<<"\n";
     std::cout<<"p1 "<<x[1]<<" "<<y[1]<<" "<<z[1]<<"\n";
     std::cout<<"p2 "<<x[2]<<" "<<y[2]<<" "<<z[2]<<"\n";
     std::cout<<"g= "<<g[0]<<" "<<g[1]<<" "<<g[2]<<"\n";
     std::cout<<"e1: "<<e1x<<" "<<e1y<<" "<<e1z<<"\n";
     exit(0);
   }
   real abse2=sqrt(e2x*e2x+e2y*e2y+e2z*e2z);
   e2x /=abse2;
   e2y /=abse2;
   e2z /=abse2;

   // check direction by using n . (e1 x e2)
   real dir= nx*(e1y*e2z-e1z*e2y) + ny*(e1z*e2x-e1x*e2z) + nz*(e1x*e2y-e1y*e2x);
   if(dir>0.0){
     e1x=-e1x;
     e1y=-e1y;
     e1z=-e1z;
   }

#if 0
   // orthogonal check
   if( fabs(nx*e1x+ny*e1y+nz*e1z)
      +fabs(nx*e2x+ny*e2y+nz*e2z)
      +fabs(e1x*e2x+e1y*e2y+e1z*e2z)>=1.0e-12){
     std::cout<<"arrange_poly:orthogonal check "<<nx*e1x+ny*e1y+nz*e1z<<" "
              <<nx*e2x+ny*e2y+nz*e2z<<" "
              <<e1x*e2x+e1y*e2y+e1z*e2z<<"\n";
   }
#endif   

   real alfa;
   for(int n=0; n<nn; n++) {
     real xtmp=x[n]-g[0];
     real ytmp=y[n]-g[1];
     real ztmp=z[n]-g[2];
     real e1=e1x*xtmp + e1y*ytmp + e1z*ztmp;
     real e2=e2x*xtmp + e2y*ytmp + e2z*ztmp;
     alfa = atan2(e1,e2);
     AngleNode an( alfa, x[n], y[n], z[n]);
     nodes.push_back(an);
   }

   /* sort the array "nodes" */
   stable_sort(nodes.begin(), nodes.end());

   /* place the sorted vales from "nodes" back to coordinates */
   for(int n=0; n<nn; n++) {
     x[n] = nodes[n].x;
     y[n] = nodes[n].y;
     z[n] = nodes[n].z;
   }
}
