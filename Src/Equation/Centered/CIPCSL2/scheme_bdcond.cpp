#include "cipcsl2.h"
using namespace std;

/******************************************************************************/
void Scheme::bdcond_i(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for a line-i variable.
*         scalar_exchange(_all) should take account of periodic condition.
*           1st: sigx.bdcond_i();
*           2nd: sigx.exchange_all();
******************************************************************************/
  Formula form;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::imin() && d != Dir::imax()) {

      /*------------------------+ 
      |  dirichlet (and inlet)  |
      +------------------------*/
      if( sca.bc().type(b) == BndType::dirichlet() ||
          sca.bc().type(b) == BndType::inlet() ) {

        int jof=1, kof=1;
        if(d == Dir::jmin() || d == Dir::jmax()) jof=0;
        if(d == Dir::kmin() || d == Dir::kmax()) kof=0;

        int jof2=0, kof2=0;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;

        /* formula is defined */
        if( sca.bc().formula(b) ) {
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()    ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            std::stringstream x, y, z, fs;
            x << "x=" << sca.xc(i); form.evaluate(x);
            y << "y=" << sca.yn(j); form.evaluate(y);
            z << "z=" << sca.zn(k); form.evaluate(z);
            fs << sca.bc().formula(b);
            real dx=sca.dxc(i);
            sigx[i][j+jof2][k+kof2] = form.evaluate(fs)*dx;
          }}}
        } else {
        /* formula is not defined */
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()    ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            real dx=sca.dxc(i);
            sigx[i][j+jof2][k+kof2] = sca.bc().value(b)*dx;
          }}}
        }
      }
      /*---------+ 
      |  others  |
      +---------*/
      if( sca.bc().type(b) == BndType::neumann()
        ||sca.bc().type(b) == BndType::symmetry()
        ||sca.bc().type(b) == BndType::wall()
        ||sca.bc().type(b) == BndType::outlet() ) {

        int jof2=0, kof2=0;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;

        if(d == Dir::jmin() || d == Dir::jmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+1; k++){
            real dx=sca.dxc(i);
            //sigx[i][j+jof2][k+kof2] = dx*0.5*(sca[i][j][k]+sca[i][j][k-1]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j][k-1]);
            stmp = min(1.0,max(0.0,stmp));
            sigx[i][j+jof2][k+kof2] = dx*stmp;
          }}}
        } else if(d == Dir::kmin() || d == Dir::kmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+1; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dx=sca.dxc(i);
            //sigx[i][j+jof2][k+kof2] = dx*0.5*(sca[i][j][k]+sca[i][j-1][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j-1][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigx[i][j+jof2][k+kof2] = dx*stmp;
          }}}
        }
      }
    }
  }

  // insert has the top priority. To overwrite the edge, call last.
  for( int b=0; b<sca.bc().count(); b++ ) {
    if( sca.bc().type_decomp(b) ) continue;
    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::imin() && d != Dir::imax()) {
      if( sca.bc().type(b) == BndType::insert() ) {

        int jof2=0, kof2=0;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;

        if(d == Dir::jmin() || d == Dir::jmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+1; k++){
            real dx=sca.dxc(i);
            //sigx[i][j+jof2][k+kof2] = dx*0.5*(sca[i][j][k]+sca[i][j][k-1]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j][k-1]);
            stmp = min(1.0,max(0.0,stmp));
            sigx[i][j+jof2][k+kof2] = dx*stmp;
          }}}
        } else if(d == Dir::kmin() || d == Dir::kmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+1; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dx=sca.dxc(i);
            //sigx[i][j+jof2][k+kof2] = dx*0.5*(sca[i][j][k]+sca[i][j-1][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j-1][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigx[i][j+jof2][k+kof2] = dx*stmp;
          }}}
        }
      }
    }
  }
}

/******************************************************************************/
void Scheme::bdcond_j(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for a line-j variable.
*         scalar_exchange(_all) should take account of periodic condition.   
*           1st: sigy.bdcond_j();                                            
*           2nd: sigy.exchange_all();                                        
******************************************************************************/
  Formula form;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::jmin() && d != Dir::jmax()) {

      /*------------------------+ 
      |  dirichlet (and inlet)  |
      +------------------------*/
      if( sca.bc().type(b) == BndType::dirichlet() ||
          sca.bc().type(b) == BndType::inlet() ) {

        int iof=1, kof=1;
        if(d == Dir::imin() || d == Dir::imax()) iof=0;
        if(d == Dir::kmin() || d == Dir::kmax()) kof=0;

        int iof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::kmin()) kof2=1;

        /* formula is defined */
        if( sca.bc().formula(b) ) {
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()    ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            std::stringstream x, y, z, fs;
            x << "x=" << sca.xc(i); form.evaluate(x);
            y << "y=" << sca.yn(j); form.evaluate(y);
            z << "z=" << sca.zn(k); form.evaluate(z);
            fs << sca.bc().formula(b);
            real dy=sca.dyc(j);
            sigy[i+iof2][j][k+kof2] = form.evaluate(fs)*dy;
          }}}
        } else {
        /* formula is not defined */
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()    ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            real dy=sca.dyc(j);
            sigy[i+iof2][j][k+kof2] = sca.bc().value(b)*dy;
          }}}
        }
      }
      /*---------+ 
      |  others  |
      +---------*/
      if( sca.bc().type(b) == BndType::neumann()
        ||sca.bc().type(b) == BndType::symmetry()
        ||sca.bc().type(b) == BndType::wall()
        ||sca.bc().type(b) == BndType::outlet() ) {

        int iof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::kmin()) kof2=1;

        if(d == Dir::imin() || d == Dir::imax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+1; k++){
            real dy=sca.dyc(j);
            //sigy[i+iof2][j][k+kof2] = dy*0.5*(sca[i][j][k]+sca[i][j][k-1]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j][k-1]);
            stmp = min(1.0,max(0.0,stmp));
            sigy[i+iof2][j][k+kof2] = dy*stmp;
          }}}
        } else if(d == Dir::kmin() || d == Dir::kmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+1; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dy=sca.dyc(j);
            //sigy[i+iof2][j][k+kof2] = dy*0.5*(sca[i][j][k]+sca[i-1][j][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i-1][j][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigy[i+iof2][j][k+kof2] = dy*stmp;
          }}}
        }
      }
    }
  }

  // insert has the top priority. To overwrite the edge, call last.
  for( int b=0; b<sca.bc().count(); b++ ) {
    if( sca.bc().type_decomp(b) ) continue;
    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::jmin() && d != Dir::jmax()) {
      if( sca.bc().type(b) == BndType::insert() ) {
        int iof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::kmin()) kof2=1;

        if(d == Dir::imin() || d == Dir::imax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+1; k++){
            real dy=sca.dyc(j);
            //sigy[i+iof2][j][k+kof2] = dy*0.5*(sca[i][j][k]+sca[i][j][k-1]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j][k-1]);
            stmp = min(1.0,max(0.0,stmp));
            sigy[i+iof2][j][k+kof2] = dy*stmp;
          }}}
        } else if(d == Dir::kmin() || d == Dir::kmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+1; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dy=sca.dyc(j);
            //sigy[i+iof2][j][k+kof2] = dy*0.5*(sca[i][j][k]+sca[i-1][j][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i-1][j][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigy[i+iof2][j][k+kof2] = dy*stmp;
          }}}
        }
      }
    }
  }
}
/******************************************************************************/
void Scheme::bdcond_k(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for a line-k variable.
*         scalar_exchange(_all) should take account of periodic condition.   
*           1st: sigz.bdcond_k();                                            
*           2nd: sigz.exchange_all();                                        
******************************************************************************/
  Formula form;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::kmin() && d != Dir::kmax()) {

      /*------------------------+ 
      |  dirichlet (and inlet)  |
      +------------------------*/
      if( sca.bc().type(b) == BndType::dirichlet() ||
          sca.bc().type(b) == BndType::inlet() ) {

        int iof=1, jof=1;
        if(d == Dir::imin() || d == Dir::imax()) iof=0;
        if(d == Dir::jmin() || d == Dir::jmax()) jof=0;

        int iof2=0, jof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;

        /* formula is defined */
        if( sca.bc().formula(b) ) {
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()    ; k++){
            std::stringstream x, y, z, fs;
            x << "x=" << sca.xc(i); form.evaluate(x);
            y << "y=" << sca.yn(j); form.evaluate(y);
            z << "z=" << sca.zn(k); form.evaluate(z);
            fs << sca.bc().formula(b);
            real dz=sca.dzc(k);
            sigz[i+iof2][j+jof2][k] = form.evaluate(fs)*dz;
          }}}
        } else {
        /* formula is not defined */
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()    ; k++){
            real dz=sca.dzc(k);
            sigz[i+iof2][j+jof2][k] = sca.bc().value(b)*dz;
          }}}
        }
      }
      /*---------+ 
      |  others  |
      +---------*/
      if( sca.bc().type(b) == BndType::neumann()
        ||sca.bc().type(b) == BndType::symmetry()
        ||sca.bc().type(b) == BndType::wall()
        ||sca.bc().type(b) == BndType::outlet() ) {

        int iof2=0, jof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;

        if(d == Dir::imin() || d == Dir::imax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+1; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dz=sca.dzc(k);
            //sigz[i+iof2][j+jof2][k] = dz*0.5*(sca[i][j][k]+sca[i][j-1][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j-1][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigz[i+iof2][j+jof2][k] = dz*stmp;
          }}}
        } else if(d == Dir::jmin() || d == Dir::jmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+1; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dz=sca.dzc(k);
            //sigz[i+iof2][j+jof2][k] = dz*0.5*(sca[i][j][k]+sca[i-1][j][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i-1][j][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigz[i+iof2][j+jof2][k] = dz*stmp;
          }}}
        }
      }
    }
  }

  // insert has the top priority. To overwrite the edge, call last.
  for( int b=0; b<sca.bc().count(); b++ ) {
    if( sca.bc().type_decomp(b) ) continue;
    Dir d = sca.bc().direction(b);
    if(d != Dir::undefined() && d != Dir::kmin() && d != Dir::kmax()) {
      if( sca.bc().type(b) == BndType::insert() ) {

        int iof2=0, jof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;

        if(d == Dir::imin() || d == Dir::imax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()  ; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+1; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dz=sca.dzc(k);
            //sigz[i+iof2][j+jof2][k] = dz*0.5*(sca[i][j][k]+sca[i][j-1][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i][j-1][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigz[i+iof2][j+jof2][k] = dz*stmp;
          }}}
        } else if(d == Dir::jmin() || d == Dir::jmax()){
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+1; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()  ; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()  ; k++){
            real dz=sca.dzc(k);
            //sigz[i+iof2][j+jof2][k] = dz*0.5*(sca[i][j][k]+sca[i-1][j][k]);
            real stmp = 0.5*(sca[i][j][k]+sca[i-1][j][k]);
            stmp = min(1.0,max(0.0,stmp));
            sigz[i+iof2][j+jof2][k] = dz*stmp;
          }}}
        }
      }
    }
  }
}

/******************************************************************************/
void Scheme::bdcond_f(const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for a nodal variable.
*         scalar_exchange(_all) should take account of periodic condition.   
*           1st: scheme_f.bdcond_i();
*           2nd: scheme_f.exchange_all();
******************************************************************************/

  Formula form;

  for( int b=0; b<sca.bc().count(); b++ ) {

    if( sca.bc().type_decomp(b) ) continue;

    /*------------------------+ 
    |  dirichlet (and inlet)  |
    +------------------------*/
    if( sca.bc().type(b) == BndType::dirichlet() ||
        sca.bc().type(b) == BndType::inlet() ) {

      Dir d = sca.bc().direction(b);
      if(d != Dir::undefined()) {
        int iof=1, jof=1, kof=1;
        if(d == Dir::imin() || d == Dir::imax()) iof=0;
        if(d == Dir::jmin() || d == Dir::jmax()) jof=0;
        if(d == Dir::kmin() || d == Dir::kmax()) kof=0;

        int iof2=0, jof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;


        /* formula is defined */
        if( sca.bc().formula(b) ) {
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            std::stringstream x, y, z, fs;
            x << "x=" << sca.xn(i); form.evaluate(x);
            y << "y=" << sca.yn(j); form.evaluate(y);
            z << "z=" << sca.zn(k); form.evaluate(z);
            fs << sca.bc().formula(b);
            f[i][j][k] = form.evaluate(fs);
          }}}
        }
        /* formula is not defined */
        else {
          for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
          for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
          for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
            f[i+iof2][j+jof2][k+kof2] = sca.bc().value(b);
          }}}
        }
      }
    }

    /*---------+ 
    |  others  |
    +---------*/
    if( sca.bc().type(b) == BndType::neumann()
      ||sca.bc().type(b) == BndType::symmetry()
      ||sca.bc().type(b) == BndType::wall()
      ||sca.bc().type(b) == BndType::outlet() ) {

      Dir d = sca.bc().direction(b);
      if(d != Dir::undefined()) {
        int iof=1, jof=1, kof=1;
        if(d == Dir::imin() || d == Dir::imax()) iof=0;
        if(d == Dir::jmin() || d == Dir::jmax()) jof=0;
        if(d == Dir::kmin() || d == Dir::kmax()) kof=0;

        int iof2=0, jof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;

        for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
        for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
        for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
          real ftmp =
                       jof*kof*0.25*(sca[i][j][k    ]+sca[i][j-jof][k    ]
                                    +sca[i][j][k-kof]+sca[i][j-jof][k-kof])
                      +iof*kof*0.25*(sca[i][j][k    ]+sca[i-iof][j][k    ]
                                    +sca[i][j][k-kof]+sca[i-iof][j][k-kof])
                      +iof*jof*0.25*(sca[i][j    ][k]+sca[i-iof][j    ][k]
                                    +sca[i][j-jof][k]+sca[i-iof][j-jof][k]);
          f[i+iof2][j+jof2][k+kof2] = min(1.0,max(0.0,ftmp));
        }}}
      }
    }
  }

  // insert has the top priority. To overwrite the edge, call last.
  for( int b=0; b<sca.bc().count(); b++ ) {
    if( sca.bc().type_decomp(b) ) continue;
    /*---------+ 
    |  insert  |
    +---------*/
    if( sca.bc().type(b) == BndType::insert() ) {

      Dir d = sca.bc().direction(b);
      if(d != Dir::undefined()) {
        int iof=1, jof=1, kof=1;
        if(d == Dir::imin() || d == Dir::imax()) iof=0;
        if(d == Dir::jmin() || d == Dir::jmax()) jof=0;
        if(d == Dir::kmin() || d == Dir::kmax()) kof=0;

        int iof2=0, jof2=0, kof2=0;
        if(d == Dir::imin()) iof2=1;
        if(d == Dir::jmin()) jof2=1;
        if(d == Dir::kmin()) kof2=1;

        for(int i=sca.bc().at(b).si(); i<=sca.bc().at(b).ei()+iof; i++){
        for(int j=sca.bc().at(b).sj(); j<=sca.bc().at(b).ej()+jof; j++){
        for(int k=sca.bc().at(b).sk(); k<=sca.bc().at(b).ek()+kof; k++){
          real ftmp =
                       jof*kof*0.25*(sca[i][j][k    ]+sca[i][j-jof][k    ]
                                    +sca[i][j][k-kof]+sca[i][j-jof][k-kof])
                      +iof*kof*0.25*(sca[i][j][k    ]+sca[i-iof][j][k    ]
                                    +sca[i][j][k-kof]+sca[i-iof][j][k-kof])
                      +iof*jof*0.25*(sca[i][j    ][k]+sca[i-iof][j    ][k]
                                    +sca[i][j-jof][k]+sca[i-iof][j-jof][k]);
          f[i+iof2][j+jof2][k+kof2] = min(1.0,max(0.0,ftmp));
        }}}
      }
    }
  }
}
