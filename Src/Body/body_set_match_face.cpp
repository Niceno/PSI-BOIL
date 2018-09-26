#include "body.h"

/******************************************************************************/
void Body::set_match_face( Cutcell & ccell, int index, int ic, int jc, int kc
                        ,const real x0, const real y0, const real z0
                        ,const real dx, const real dy, const real dz
                        ,int face_match[][2], real face_match_norm[][2][3]
                        ,int & nface_match) {

  real xp[4];
  real yp[4];
  real zp[4];


    std::cout<<"cut_cell: nface_match=b "<<ic<<" "<<jc<<" "<<nface_match<<"\n";
    if(face_match[0][0]==1&&face_match_norm[0][0][0]<0.0) return;
    if(face_match[0][1]==1&&face_match_norm[0][1][0]>0.0) return;
    if(face_match[1][0]==1&&face_match_norm[1][0][1]<0.0) return;
    if(face_match[1][1]==1&&face_match_norm[1][1][1]>0.0) return;
    if(face_match[2][0]==1&&face_match_norm[2][0][2]<0.0) return;
    if(face_match[2][1]==1&&face_match_norm[2][1][2]>0.0) return;
      // initialize as 1.0
      ccell->fS(1.0,W_E,0);
      ccell->fS(1.0,W_E,1);
      ccell->fS(1.0,S_N,0);
      ccell->fS(1.0,S_N,1);
      ccell->fS(1.0,B_T,0);
      ccell->fS(1.0,B_T,1);
    std::cout<<"cut_cell: nface_match=a "<<ic<<" "<<jc<<" "<<nface_match<<"\n";
    int nn=4;
    /* i-min */
    if(face_match[0][0]==1){
      xp[0]=x0;  yp[0]=y0;     zp[0]=z0;
      xp[1]=x0;  yp[0]=y0+dy;  zp[0]=z0;
      xp[2]=x0;  yp[2]=y0+dy;  zp[2]=z0+dz;
      xp[3]=x0;  yp[0]=y0;     zp[0]=z0+dz;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,W_E,0);
    }
    if(face_match[0][1]==1){
      ccell->fS(0.0,W_E,1);
    }
    if(face_match[1][0]==1){
      ccell->fS(0.0,S_N,0);
    }
    if(face_match[1][1]==1){
      ccell->fS(0.0,S_N,1);
    }
    if(face_match[2][0]==1){
      ccell->fS(0.0,B_T,0);
    }
    if(face_match[2][1]==1){
      ccell->fS(0.0,B_T,1);
    }
    ccell->fV(1.0);
    return;

}

/*-----------------------------------------------------------------------------+
 '$Id: body_set_match_face.cpp,v 1.1 2014/02/03 14:12:34 sato Exp $'/
+-----------------------------------------------------------------------------*/
