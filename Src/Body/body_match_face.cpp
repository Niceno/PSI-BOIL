#include "body.h"

/******************************************************************************/
void Body::match_face( int index, int ic, int jc, int kc
                        ,real x[], real y[], real z[]
                        ,int face_match[][2], real face_match_norm[][2][3]
                        ,int & nface_match) {
#if 0
  if(ic==20&&jc==19&&kc==5){
    std::cout<<"face_match\n";
  }
#endif
  /*----------------------------------------+ 
  |  do not browse through all polygons,    |
  |  only throught those stored in polytag  |
  +----------------------------------------*/
  int nni[2]={0,0};
  int nnj[2]={0,0};
  int nnk[2]={0,0};
  for(int p=0; p<polytags[index].size(); p++) {
    const int c = polytags[index][p];

    // x-direction
    if( fabs(polys[c].n(0))> 0.999 ) {
      for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
          for(int k=0; k<2; k++){
            real xi = polys[c].a()*x[i]
                    + polys[c].b()*y[j]
                    + polys[c].c()*z[k]
                    + polys[c].d();
            if(fabs(xi)<tol){
               real A0 = area(y[j], z[k],
                         polys[c].yn(1), polys[c].zn(1),
                         polys[c].yn(2), polys[c].zn(2));

               real A1 = area(polys[c].yn(0), polys[c].zn(0),
                         y[j], z[k],
                         polys[c].yn(2), polys[c].zn(2));

               real A2 = area(polys[c].yn(0), polys[c].zn(0),
                         polys[c].yn(1), polys[c].zn(1),
                         y[j], z[k]);

               bool inside = false;
               if( polys[c].n(0)>-boil::pico && A0>=-boil::pico
                          && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
               if( polys[c].n(0)< boil::pico && A0<= boil::pico
                          && A1<= boil::pico && A2<= boil::pico ) inside=true;
               if(inside){
                 nni[i]++;
                 face_match_norm[0][i][0]+=polys[c].n(0);
                 face_match_norm[0][i][1]+=polys[c].n(1);
                 face_match_norm[0][i][2]+=polys[c].n(2);
               }
            }
          }
        }
      }
    }

    // y-direction
    if( fabs(polys[c].n(1))> 0.999 ) {
      for(int j=0; j<2; j++){
        for(int i=0; i<2; i++){
          for(int k=0; k<2; k++){
            real xi = polys[c].a()*x[i]
                    + polys[c].b()*y[j]
                    + polys[c].c()*z[k]
                    + polys[c].d();
            if(fabs(xi)<tol){
              real A0 = area(z[k], x[i],
                             polys[c].zn(1), polys[c].xn(1),
                             polys[c].zn(2), polys[c].xn(2));

              real A1 = area(polys[c].zn(0), polys[c].xn(0),
                             z[k], x[i],
                             polys[c].zn(2), polys[c].xn(2));

              real A2 = area(polys[c].zn(0), polys[c].xn(0),
                             polys[c].zn(1), polys[c].xn(1),
                             z[k], x[i]);

              bool inside = false;
              if( polys[c].n(1)>-boil::pico && A0>=-boil::pico
                         && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
              if( polys[c].n(1)< boil::pico && A0<= boil::pico
                         && A1<= boil::pico && A2<= boil::pico ) inside=true;

               if(inside){
                 nnj[j]++;
                 face_match_norm[1][j][0]+=polys[c].n(0);
                 face_match_norm[1][j][1]+=polys[c].n(1);
                 face_match_norm[1][j][2]+=polys[c].n(2);
               }
            }
          }
        }
      }
    }

    // z-direction
    if( fabs(polys[c].n(2))> 0.999 ) {
      for(int k=0; k<2; k++){
        for(int i=0; i<2; i++){
          for(int j=0; j<2; j++){
            real xi = polys[c].a()*x[i]
                    + polys[c].b()*y[j]
                    + polys[c].c()*z[k]
                    + polys[c].d();
            if(fabs(xi)<tol){
              real A0 = area(x[i], y[j],
                             polys[c].xn(1), polys[c].yn(1),
                             polys[c].xn(2), polys[c].yn(2));

              real A1 = area(polys[c].xn(0), polys[c].yn(0),
                             x[i], y[j],
                             polys[c].xn(2), polys[c].yn(2));
      
              real A2 = area(polys[c].xn(0), polys[c].yn(0),
                             polys[c].xn(1), polys[c].yn(1),
                             x[i], y[j]);

              bool inside = false;
              if( polys[c].n(2)>-boil::pico && A0>=-boil::pico
                         && A1>=-boil::pico && A2>=-boil::pico ) inside=true;
              if( polys[c].n(2)< boil::pico && A0<= boil::pico
                         && A1<= boil::pico && A2<= boil::pico ) inside=true;

              if(inside){
                nnk[k]++;
                 face_match_norm[2][k][0]+=polys[c].n(0);
                 face_match_norm[2][k][1]+=polys[c].n(1);
                 face_match_norm[2][k][2]+=polys[c].n(2);
              }
            }
          }
        }
      }
    }
  }

  for(int i=0; i<2; i++){
    if(nni[i]>=4){
      face_match[0][i]=1;
      face_match_norm[0][i][0]/=real(nni[i]);
      face_match_norm[0][i][1]/=real(nni[i]);
      face_match_norm[0][i][2]/=real(nni[i]);
      real mag=sqrt( pow(face_match_norm[0][i][0],2.0)
                   + pow(face_match_norm[0][i][1],2.0)
                   + pow(face_match_norm[0][i][2],2.0));
      face_match_norm[0][i][0]/=mag;
      face_match_norm[0][i][1]/=mag;
      face_match_norm[0][i][2]/=mag;
    }
  }
  for(int j=0; j<2; j++){
    if(nnj[j]>=4){
      face_match[1][j]=1;
      face_match_norm[1][j][0]/=real(nnj[j]);
      face_match_norm[1][j][1]/=real(nnj[j]);
      face_match_norm[1][j][2]/=real(nnj[j]);
      real mag=sqrt( pow(face_match_norm[1][j][0],2.0)
                   + pow(face_match_norm[1][j][1],2.0)
                   + pow(face_match_norm[1][j][2],2.0));
      face_match_norm[1][j][0]/=mag;
      face_match_norm[1][j][1]/=mag;
      face_match_norm[1][j][2]/=mag;
    }
  }
  for(int k=0; k<2; k++){
    if(nnk[k]>=4){
      face_match[2][k]=1;
      face_match_norm[2][k][0]/=real(nnk[k]);
      face_match_norm[2][k][1]/=real(nnk[k]);
      face_match_norm[2][k][2]/=real(nnk[k]);
      real mag=sqrt( pow(face_match_norm[2][k][0],2.0)
                   + pow(face_match_norm[2][k][1],2.0)
                   + pow(face_match_norm[2][k][2],2.0));
      face_match_norm[2][k][0]/=mag;
      face_match_norm[2][k][1]/=mag;
      face_match_norm[2][k][2]/=mag;
    }
  }

  for(int i=0; i<3; i++){
    for(int j=0; j<2; j++){
      if(face_match[i][j]==1){
        nface_match++;
#if 0
  if(ic==20&&jc==19&&kc==5){
    std::cout<<"face_match: "<<i<<" "<<j<<"\n";
  }
#endif
      }
    }
  }

}

/*-----------------------------------------------------------------------------+
 '$Id: body_match_face.cpp,v 1.1 2014/02/03 14:12:33 sato Exp $'/
+-----------------------------------------------------------------------------*/
