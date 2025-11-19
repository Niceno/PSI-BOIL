#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"

/******************************************************************************/
real cal_area(real vec1[], real vec2[]) {
  real ox = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  real oy = vec1[2]*vec2[0] - vec1[0]*vec2[2];
  real oz = vec1[0]*vec2[1] - vec1[1]*vec2[0];
  return 0.5*sqrt(ox*ox+oy*oy+oz*oz);
}

/******************************************************************************/
real face_area(real xp[], real yp[], real zp[], const int n) {

  assert( n>=3 );
  real vec1[3],vec2[3];
  real area=0.0;

  for(int i=0; i<=n-3; i++){
    vec1[0]=xp[i+1]-xp[0];
    vec1[1]=yp[i+1]-yp[0];
    vec1[2]=zp[i+1]-zp[0];
    vec2[0]=xp[i+2]-xp[0];
    vec2[1]=yp[i+2]-yp[0];
    vec2[2]=zp[i+2]-zp[0];
    area += cal_area(vec1, vec2);
  }
  return area;
}

/******************************************************************************/
void Body::fluid_area(CutCell *ccell, int n_in_fluid[][2][2],
                bool edge_cut[][2][2],
                real x_[][2], real y_[][2], real z_[][2],
                real x[], real y[], real z[],
                real dx,  real dy,  real dz) {

  const int W_E = 0;
  const int S_N = 1;
  const int B_T = 2;
  const real tolv=1e-8;

  int n;
  real xp[8];
  real yp[8];
  real zp[8];

  /*----------------------------------+
  |  determine if faces are in fluid  | -> uses information from nodes
  +----------------------------------*/
  int f_in_fluid[3][2];
        
  /* w - e */
  for(int i=0; i<2; i++)
    f_in_fluid[W_E][i] = int(ccell->fP(i,0,0) + ccell->fP(i,1,0) +
                             ccell->fP(i,0,1) + ccell->fP(i,1,1));
          
  /* s - n */
  for(int j=0; j<2; j++)
    f_in_fluid[S_N][j] = int(ccell->fP(0,j,0) + ccell->fP(1,j,0) +
                             ccell->fP(0,j,1) + ccell->fP(1,j,1));
        
  /* b - t */
  for(int k=0; k<2; k++)
    f_in_fluid[B_T][k] = int(ccell->fP(0,0,k) + ccell->fP(1,0,k) +
                             ccell->fP(0,1,k) + ccell->fP(1,1,k));

  /*--------------------+
  |  face area in fluid |
  +--------------------*/
  for(int off=0; off<2; off++){
    for(int d=0; d<3; d++) {
      switch(d) {
        case W_E: if( f_in_fluid[W_E][off]==4 ){
                  ccell->fS(1.0,W_E,off);
                } else if (f_in_fluid[W_E][off]==0){
                  ccell->fS(0.0,W_E,off);
                } else {
                  int iflag = n_in_fluid[off][0][0] + n_in_fluid[off][0][1]
                            + n_in_fluid[off][1][0] + n_in_fluid[off][1][1];
                  if(iflag==8){         //2+2+2+2
                    ccell->fS(1.0,W_E,off);
                  } else if(iflag==6){         //2+2+1+1
                    ccell->fS(1.0,W_E,off);
                  } else if(iflag==5){  //2+2+1+0
                    ccell->fS(0.5-tolv,W_E,off);
                  } else {
                    n=0;
                    if(n_in_fluid[off][0][0]>=1){
                      xp[n]=x[off]; yp[n]=y[0]; zp[n]=z[0]; n++;
                    }
                    if(edge_cut[1][off][0]){
                      xp[n]=x[off]; yp[n]=y_[off][0]; zp[n]=z[0]; n++;
                    }
                    if(n_in_fluid[off][1][0]>=1){
                      xp[n]=x[off]; yp[n]=y[1]; zp[n]=z[0]; n++;
                    }
                    if(edge_cut[2][off][1]){
                      xp[n]=x[off]; yp[n]=y[1]; zp[n]=z_[off][1]; n++;
                    }
                    if(n_in_fluid[off][1][1]>=1){
                      xp[n]=x[off]; yp[n]=y[1]; zp[n]=z[1]; n++;
                    }
                    if(edge_cut[1][off][1]){
                      xp[n]=x[off]; yp[n]=y_[off][1]; zp[n]=z[1]; n++;
                    }
                    if(n_in_fluid[off][0][1]>=1){
                      xp[n]=x[off]; yp[n]=y[0]; zp[n]=z[1]; n++;
                    }
                    if(edge_cut[2][off][0]){
                      xp[n]=x[off]; yp[n]=y[0]; zp[n]=z_[off][0]; n++;
                    }
                    assert(n>=3);
                    ccell->fS(face_area(xp, yp, zp, n)/(dy*dz),W_E,off);
                  }
                }
                break;
        case S_N: if( f_in_fluid[S_N][off]==4 ){
                  ccell->fS(1.0,S_N,off);
                } else if ( f_in_fluid[S_N][off]==0 ){
                  ccell->fS(0.0,S_N,off);
                } else { 
                  int iflag = n_in_fluid[0][off][0] + n_in_fluid[0][off][1]
                            + n_in_fluid[1][off][0] + n_in_fluid[1][off][1];
                  if(iflag==8){         //2+2+2+2
                    ccell->fS(1.0,S_N,off); 
                  } else if(iflag==6){         //2+2+1+1
                    ccell->fS(1.0,S_N,off); 
                  } else if(iflag==5){  //2+2+1+0 
                    ccell->fS(0.5-tolv,S_N,off);
                  } else { 
                    n=0;
                    if(n_in_fluid[0][off][0]>=1){
                      xp[n]=x[0]; yp[n]=y[off]; zp[n]=z[0]; n++;
                    }
                    if(edge_cut[0][off][0]){ 
                      xp[n]=x_[off][0]; yp[n]=y[off]; zp[n]=z[0]; n++;
                    }
                    if(n_in_fluid[1][off][0]>=1){
                      xp[n]=x[1]; yp[n]=y[off]; zp[n]=z[0]; n++;
                    }
                    if(edge_cut[2][1][off]){
                      xp[n]=x[1]; yp[n]=y[off]; zp[n]=z_[1][off]; n++;
                    }
                    if(n_in_fluid[1][off][1]>=1){
                      xp[n]=x[1]; yp[n]=y[off]; zp[n]=z[1]; n++;
                    }
                    if(edge_cut[0][off][1]){
                      xp[n]=x_[off][1]; yp[n]=y[off]; zp[n]=z[1]; n++;
                    }
                    if(n_in_fluid[0][off][1]>=1){
                      xp[n]=x[0]; yp[n]=y[off]; zp[n]=z[1]; n++;
                    }
                    if(edge_cut[2][0][off]){
                      xp[n]=x[0]; yp[n]=y[off]; zp[n]=z_[0][off]; n++;
                    }
                    assert(n>=3);
                    ccell->fS(face_area(xp, yp, zp, n)/(dx*dz),S_N,off);
                  }
                }
                break;
        case B_T: if( f_in_fluid[B_T][off]==4 ){
                  ccell->fS(1.0,B_T,off);
                } else if( f_in_fluid[B_T][off]==0 ){
                  ccell->fS(0.0,B_T,off);
                } else {
                  int iflag = n_in_fluid[0][0][off] + n_in_fluid[1][0][off]
                            + n_in_fluid[0][1][off] + n_in_fluid[1][1][off];
                  if(iflag==8){         //2+2+2+2
                    ccell->fS(1.0,B_T,off);
                  } else if(iflag==6){         //2+2+1+1
                    ccell->fS(1.0,B_T,off);
                  } else if(iflag==5){  //2+2+1+0
                    ccell->fS(0.5-tolv,B_T,off);
                  } else {
                    n=0;
                    if(n_in_fluid[0][0][off]>=1){
                      xp[n]=x[0]; yp[n]=y[0]; zp[n]=z[off]; n++;
                    }
                    if(edge_cut[1][0][off]){
                      xp[n]=x[0]; yp[n]=y_[0][off]; zp[n]=z[off]; n++;
                    }
                    if(n_in_fluid[0][1][off]>=1){
                      xp[n]=x[0]; yp[n]=y[1]; zp[n]=z[off]; n++;
                    }
                    if(edge_cut[0][1][off]){
                      xp[n]=x_[1][off]; yp[n]=y[1]; zp[n]=z[off]; n++;
                    }
                    if(n_in_fluid[1][1][off]>=1){
                      xp[n]=x[1]; yp[n]=y[1]; zp[n]=z[off]; n++;
                    }
                    if(edge_cut[1][1][off]){
                      xp[n]=x[1]; yp[n]=y_[1][off]; zp[n]=z[off]; n++;
                    }
                    if(n_in_fluid[1][0][off]>=1){
                      xp[n]=x[1]; yp[n]=y[0]; zp[n]=z[off]; n++;
                    }
                    if(edge_cut[0][0][off]){
                      xp[n]=x_[0][off]; yp[n]=y[0]; zp[n]=z[off]; n++;
                    }
                    assert(n>=3);
                    ccell->fS(face_area(xp, yp, zp, n)/(dx*dy),B_T,off);
                  }
                }
                break;

      }
    }
  }
  return ;
}
