#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_norm_cc(const Scalar & val) {
/***************************************************************************//**
*  \brief normal vector for cells adjacent to a wall or an immersed boundary
*******************************************************************************/

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      /*-------+
      |  Wall  |
      +-------*/

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1;
            norm_cc_imin(val,ii,j,k);
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;
            norm_cc_imax(val,ii,j,k);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;
            norm_cc_jmin(val,i,jj,k);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1;
            norm_cc_jmax(val,i,jj,k);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1;
            norm_cc_kmin(val,i,j,kk);
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1;
            norm_cc_kmax(val,i,j,kk);
          }
        }
      }
    }

  } /* bcs */

  /***************+
  | immersed body |
  +***************/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    // cell[i][j][k] is wall adjacent cells in fluid domain
    dom->ibody().ijk(cc,&i,&j,&k);

    // west is in solid domain
    if (dom->ibody().off(i-1,j,k)) {
      norm_cc_imin(val,i,j,k);
    }

    // east
    if (dom->ibody().off(i+1,j,k)) {
      norm_cc_imax(val,i,j,k);
    }

    // south
    if (dom->ibody().off(i,j-1,k)) {
      norm_cc_kmin(val,i,j,k);
    }

    // north
    if (dom->ibody().off(i,j+1,k)) {
      norm_cc_jmax(val,i,j,k);
    }

    // bottom
    if (dom->ibody().off(i,j,k-1)) {
      norm_cc_kmin(val,i,j,k);
    }

    // top
    if (dom->ibody().off(i,j,k+1)) {
      norm_cc_kmax(val,i,j,k);
    }
  }

  /*****************+
  | corner of walls |
  +*****************/
  /* line i-min & j-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=si();
    int j=sj();
    for_k(k){
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]));
      nzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));

      nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-min & j-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    for_k(k) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      nzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-min & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int k=sk();
    for_j(j) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      nyX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]));

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));

      nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-min & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int k=ek();
    for_j(j) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      nyX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i  ][j][k-1]+val[i  ][j][k]));
      nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-max & j-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    for_k(k) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j  ][k]+val[i-1][j  ][k]));
      nzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-max & j-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    for_k(k) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      nzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-max & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int k=sk();
    for_j(j) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      nyX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      nzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k  ]+val[i-1][j][k  ]));

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));

      nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line i-max & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int k=ek();
    for_j(j) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      nyX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      nzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line j-min & k-min */
  if(val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=sj();
    int k=sk();
    for_i(i) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
      nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));

      nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k  ]+val[i][j+1][k  ]));

      nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line j-min & k-max */
  if(val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=sj();
    int k=ek();
    for_i(i) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
      nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line j-max & k-min */
  if(val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej();
    int k=sk();
    for_i(i) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
      nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));

      nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
                  -(val[i][j-1][k  ]+val[i][j][k  ]));

      nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* line j-max & k-max */
  if(val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej();
    int k=ek();
    for_i(i) {
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
      nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
                  -(val[i][j-1][k-1]+val[i][j][k-1]));

      nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
    }
  }

  /* corner i-max & j-min & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=sk();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    nyX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    nzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-max & j-min & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=ek();
      real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
      nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      nyX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j  ][k]+val[i-1][j  ][k]));
      nzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));

      nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

      Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-min & j-min & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int j=sj();
    int k=sk();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                -(val[i+1][j  ][k]+val[i][j  ][k]));
    nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                -(val[i+1][j][k  ]+val[i][j][k  ]));

    nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
               - (val[i][j][k  ]+val[i][j+1][k  ]));

    nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-min & j-min & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int j=sj();
    int k=ek();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                -(val[i+1][j  ][k]+val[i][j  ][k]));
    nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                -(val[i+1][j][k-1]+val[i][j][k-1]));

    nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
               - (val[i][j][k-1]+val[i][j+1][k-1]));

    nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-min & j-max & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    int k=sk();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                -(val[i+1][j-1][k]+val[i][j-1][k]));
    nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                -(val[i+1][j][k  ]+val[i][j][k  ]));

    nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));

    nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-min & j-max & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    int k=ek();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                -(val[i+1][j-1][k]+val[i][j-1][k]));
    nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                -(val[i+1][j][k-1]+val[i][j][k-1]));

    nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));

    nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-max & j-min & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=sk();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    nyX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    nzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-max & j-min & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=ek();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    nyX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    nzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                -(val[i][j][k-1]+val[i-1][j][k-1]));

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                -(val[i][j][k-1]+val[i][j+1][k-1]));

    nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-max & j-max & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    int k=sk();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    nyX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                -(val[i][j-1][k]+val[i-1][j-1][k]));
    nzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));

    nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }

  /* corner i-max & j-max & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    int k=ek();
    real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
    nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    nyX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                -(val[i][j-1][k]+val[i-1][j-1][k]));
    nzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                -(val[i][j][k-1]+val[i-1][j][k-1]));

    nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));

    nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

    Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
  }
}


/******************************************************************************/
void VOF::norm_cc_imin(const Scalar & val,
                       const int i,const int j, const int k) {
  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
  nyX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]));
  nzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]));

  nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));

  nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}

/******************************************************************************/
void VOF::norm_cc_imax(const Scalar & val,
                       const int i,const int j, const int k) {
  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
  nyX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i][j-1][k]+val[i-1][j-1][k]));
  nzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i][j][k-1]+val[i-1][j][k-1]));

  nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
             - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));

  nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}

/******************************************************************************/
void VOF::norm_cc_jmin(const Scalar & val,
                       const int i,const int j, const int k) {
  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  nyX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
  nzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));

  nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
  nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j][k-1]+val[i][j+1][k-1]));

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}

/******************************************************************************/
void VOF::norm_cc_jmax(const Scalar & val,
                       const int i,const int j, const int k) {
  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  nyX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  nzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]));
  nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
              -(val[i][j-1][k-1]+val[i][j][k-1]));

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}

/******************************************************************************/
void VOF::norm_cc_kmin(const Scalar & val,
                       const int i,const int j, const int k) {
  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  nyX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  nzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));

  nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}

/******************************************************************************/
void VOF::norm_cc_kmax(const Scalar & val,
                       const int i,const int j, const int k) {

  real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;
  nxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  nyX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  nzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
              -(val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
              -(val[i-1][j][k-1]+val[i-1][j][k]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
              -(val[i][j-1][k-1]+val[i][j-1][k]));
  nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));

  Comp mcomp; 
  select_norm_cc(nx[i][j][k],ny[i][j][k],nz[i][j][k],
                 nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ, mcomp);
}
