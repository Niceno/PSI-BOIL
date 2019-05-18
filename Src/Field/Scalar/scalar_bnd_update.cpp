#include "scalar.h"

/******************************************************************************/
void Scalar::bnd_update() {
/*-----------------------------------------------------------------------------+
|  updates boundary condition values for a scalar variable.                    |
|                                                                              |
|  it doesn't work for convective boundary conditions, because the necessary   |
|  information on physical properties and transfer coeffients are missing.     |
|                                                                              |
|  scalar_exchange_all should take account of periodic condition               |
|    1st: phi.bnd_update();                                                    |
|    2nd: phi.exchange_all();                                                  |
+-----------------------------------------------------------------------------*/

  Formula F;

  for( int b=0; b<bc().count(); b++ ) {

    if( bc().type_decomp(b) ) continue;

    // set loop range
    int ist,ied,jst,jed,kst,ked;
    if(bc().at(b).si()==bc().at(b).ei()){
      if(bc().at(b).si()>boil::BW) {
        ist=ei()+1;
        ied=ei()+boil::BW;
      } else {
        ist=si()-boil::BW;
        ied=si()-1;
      }
    } else {
      if(bc().at(b).si()==si()){ist=si()-boil::BW;}
      else{ist=bc().at(b).si();}
      if(bc().at(b).ei()==ei()){ied=ei()+boil::BW;}
      else{ied=bc().at(b).ei();}
    }
    if(bc().at(b).sj()==bc().at(b).ej()){
      if(bc().at(b).sj()>boil::BW) {
        jst=ej()+1;
        jed=ej()+boil::BW;
      } else {
        jst=sj()-boil::BW;
        jed=sj()-1;
      }
    } else {
      if(bc().at(b).sj()==sj()){jst=sj()-boil::BW;}
      else{jst=bc().at(b).sj();}
      if(bc().at(b).ej()==ej()){jed=ej()+boil::BW;}
      else{jed=bc().at(b).ej();}
    }
    if(bc().at(b).sk()==bc().at(b).ek()){
      if(bc().at(b).sk()>boil::BW) {
        kst=ek()+1;
        ked=ek()+boil::BW;
      } else {
        kst=sk()-boil::BW;
        ked=sk()-1;
      }
    } else {
      if(bc().at(b).sk()==sk()){kst=sk()-boil::BW;}
      else{kst=bc().at(b).sk();}
      if(bc().at(b).ek()==ek()){ked=ek()+boil::BW;}
      else{ked=bc().at(b).ek();}
    }

    /*========================+ 
    |  dirichlet (and inlet)  |
    +========================*/
    if( bc().type(b) == BndType::dirichlet() ||
        bc().type(b) == BndType::inlet() ) {

      /* formula is defined */
      if( bc().formula(b) ) {
        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(int i=ist; i<=ied; i++)
          for(int j=jst; j<=jed; j++)
            for(int k=kst; k<=ked; k++) {
              std::stringstream x, y, z, f;
              x << "x=" << xc(i); F.evaluate(x);
              y << "y=" << yc(j); F.evaluate(y);
              z << "z=" << zc(k); F.evaluate(z);
              f << bc().formula(b);

              val[i][j][k] = F.evaluate(f);
        }
      }
      /* formula is not defined */
      else {
        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(int i=ist; i<=ied; i++)
          for(int j=jst; j<=jed; j++)
            for(int k=kst; k<=ked; k++) {
              val[i][j][k] = bc().value(b);
            } 
      }
    }

    /*==========+ 
    |  others   |
    +==========*/
    if( bc().type(b) == BndType::neumann()  ||
        bc().type(b) == BndType::wall()     ||
        bc().type(b) == BndType::outlet() ) {

      Dir d = bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      
      if( d == Dir::ibody() ) {
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) ) 
           val[i-1][j][k] = val[i][j][k];
          if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i-1][j][k];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) ) 
           val[i+1][j][k] = val[i][j][k];
          if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i+1][j][k];
  

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) ) 
           val[i][j-1][k] = val[i][j][k];
          if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j-1][k];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) ) 
           val[i][j+1][k] = val[i][j][k];
          if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j+1][k];
          
  
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) ) 
           val[i][j][k-1] = val[i][j][k];
          if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j][k-1];
          
          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) ) 
           val[i][j][k+1] = val[i][j][k];
          if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) ) 
           val[i][j][k] = val[i][j][k+1];
        }

      } else if(d != Dir::undefined()) {
        int i,j,k;
        int * iref = &i;
        int * jref = &j;
        int * kref = &k;
        int ibeg = si(), iend = ei(), 
            jbeg = sj(), jend = ej(),
            kbeg = sk(), kend = ek();

        if(d == Dir::imin()) iref = &ibeg; if(d == Dir::imax()) iref = &iend;
        if(d == Dir::jmin()) jref = &jbeg; if(d == Dir::jmax()) jref = &jend;
        if(d == Dir::kmin()) kref = &kbeg; if(d == Dir::kmax()) kref = &kend;

        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++)
              val[i][j][k]=val[*iref][*jref][*kref];
      }
    }

    /*====================+ 
    |  symmetry boundary  |
    +====================*/
    if( bc().type(b) == BndType::symmetry() ) {

      Dir d = bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/

      if( d == Dir::ibody() ) {
        for(int cc=0; cc<dom->ibody().nccells(); cc++) {
          int i,j,k;
          dom->ibody().ijk(cc,&i,&j,&k); // OPR(i); OPR(j); OPR(k);

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) )
           val[i-1][j][k] = val[i][j][k];
          if( dom->ibody().on(i-1,j,k) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i-1][j][k];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) )
           val[i+1][j][k] = val[i][j][k];
          if( dom->ibody().on(i+1,j,k) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i+1][j][k];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) )
           val[i][j-1][k] = val[i][j][k];
          if( dom->ibody().on(i,j-1,k) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i][j-1][k];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) )
           val[i][j+1][k] = val[i][j][k];
          if( dom->ibody().on(i,j+1,k) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i][j+1][k];


          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) )
           val[i][j][k-1] = val[i][j][k];
          if( dom->ibody().on(i,j,k-1) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i][j][k-1];

          if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) )
           val[i][j][k+1] = val[i][j][k];
          if( dom->ibody().on(i,j,k+1) && dom->ibody().off(i,j,k) )
           val[i][j][k] = val[i][j][k+1];
        }

      } else if(d != Dir::undefined()) {

        /* by default, the copied values are taken with the same i,j,k */
        int i,j,k;
        int * iref = &i;
        int * jref = &j;
        int * kref = &k;

        /* here, indices of domain boundary are stored,
           they will be in/decremented once before they are used for the first
           time */
        int ibeg = si()-1, iend = ei()+1,
            jbeg = sj()-1, jend = ej()+1,
            kbeg = sk()-1, kend = ek()+1;

        /* by default, for loop goes in ascending order */
        int iinc(1), jinc(1), kinc(1);
        
        /* for readability, order of for loops is decided through a flag */
        bool idir(false), jdir(false), kdir(false);

        /* if dir is imin, iref has to be shifted to ibeg. We are dealing
           with symmetry, so this index will start at si() and progress towards
           si()+boil::BW-1. At the same time, i-index should decrease from 
           ied = si()-1 to ist = si()-boil::BW = 0. As a result, the iref
           pointer is changed, the loop increment and bounds are inverted.
           The implementation of the for loop then mirrors
                         val[si()-ii][j][k] = val[si()-1+ii][j][k] 
           with ii going from 1 to boil::BW (buffer width) */
        if       (d == Dir::imin()) {
          iref = &ibeg;
          iinc = -1;
          std::swap(ist,ied);
          idir = true;
        } else if(d == Dir::imax()) {
          iref = &iend;
          idir = true;
        } else if(d == Dir::jmin()) {
          jref = &jbeg; 
          jinc = -1;
          std::swap(jst,jed);
          jdir = true;
        } else if(d == Dir::jmax()) {
          jref = &jend;
          jdir = true;
        } else if(d == Dir::kmin()) {
          kref = &kbeg; 
          kinc = -1;
          std::swap(kst,ked);
          kdir = true;
        } else if(d == Dir::kmax()) {
          kref = &kend;
          kdir = true;
        }

        /* used to be: for_vijk( bc().at(b), i,j,k ) */
        /* changes of _beg and _end only affect the direction in question */
        if(idir) {
          for(i=ist; i!=ied+iinc; i += iinc) {
            ibeg += 1;
            iend -= 1;
            for(j=jst; j!=jed+jinc; j += jinc) {
              for(k=kst; k!=ked+kinc; k += kinc) {
                val[i][j][k]=val[*iref][*jref][*kref];
              }
            }
          }
        } else if (jdir) {
          for(j=jst; j!=jed+jinc; j += jinc) {
            jbeg += 1;
            jend -= 1;
            for(i=ist; i!=ied+iinc; i += iinc) {
              for(k=kst; k!=ked+kinc; k += kinc) {
                val[i][j][k]=val[*iref][*jref][*kref];
              }
            }
          }
        } else {
          for(k=kst; k!=ked+kinc; k += kinc) {
            kbeg += 1;
            kend -= 1; 
            for(i=ist; i!=ied+iinc; i += iinc) {
              for(j=jst; j!=jed+jinc; j += jinc) {
                val[i][j][k]=val[*iref][*jref][*kref];
              }
            }
          }
        } /* which direction is looped through */
      } /* dir not undefined */
    } /* symmetry bc */

    /*======================+ 
    |  convective boundary  |
    +======================*/
    if( bc().type(b) == BndType::convective() ) {
      boil::oout<<"Obsolete boundary condition (convective)! Exiting."
                <<boil::endl;
      exit(0);
      for_vijk( bc().at(b), i,j,k ) {
        /* do nothing here, it is handled in the Centered class */
      }
    }

  }
}
