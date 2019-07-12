#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_norm_cc(const Scalar & val) {
/***************************************************************************//**
*  \brief normal vector for the adjacent cells next wall, symmetric and
*         immersed boundary
*******************************************************************************/

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::wall() ) {

      //std::cout<<"insert_bc_norm_cc: "<<boil::cart.iam()<<"\n";

      /*-------+
      |  Wall  |
      +-------*/

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;

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

#if 0
    // Symmetry can be calculated in normal way in norm_cc
    if( val.bc().type(b) == BndType::symmetry() ) {

      /*-----------+
      |  Symmetry  |
      +-----------*/

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        real nxX, nyX, nzX, nxY, nyY, nzY, nxZ, nyZ, nzZ;

        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1; // val[ii-1]=val[ii]
            nxX = copysign(1.0,+(val[ii+1][j][k]-val[ii][j][k]));
            nyX = 0.5 * ( (val[ii+1][j+1][k]+val[ii][j+1][k]+val[ii][j+1][k])
                        - (val[ii+1][j-1][k]+val[ii][j-1][k]+val[ii][j-1][k]));
            nzX = 0.5 * ( (val[ii+1][j][k+1]+val[ii][j][k+1]+val[ii][j][k+1])
                        - (val[ii+1][j][k-1]+val[ii][j][k-1]+val[ii][j][k-1]));
            normalize(nxX,nyX,nzX);

            nxY = 1.0 * ( (val[ii+1][j-1][k]+val[ii+1][j][k]+val[ii+1][j+1][k])
                        - (val[ii  ][j-1][k]+val[ii  ][j][k]+val[ii  ][j+1][k]));
            nyY = copysign(1.0,+(val[ii][j+1][k]-val[ii][j-1][k]));
            nzY = 0.5 * ( (val[ii][j-1][k+1]+val[ii][j][k+1]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j][k-1]+val[ii][j+1][k-1]));
            normalize(nxY,nyY,nzY);

            nxZ = 1.0 * ( (val[ii+1][j][k-1]+val[ii+1][j][k]+val[ii+1][j][k+1])
                        - (val[ii  ][j][k-1]+val[ii  ][j][k]+val[ii  ][j][k+1]));
            nyZ = 0.5 * ( (val[ii][j+1][k-1]+val[ii][j+1][k]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j-1][k]+val[ii][j-1][k+1]));
            nzZ = copysign(1.0,+(val[ii][j][k+1]-val[ii][j][k-1]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,ii,j,k);
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;  // val[ii+1]=val[ii]
            nxX = copysign(1.0,+(val[ii][j][k]-val[ii-1][j][k]));
            nyX = 0.5 * ( (val[ii][j+1][k]+val[ii][j+1][k]+val[ii-1][j+1][k])
                        - (val[ii][j-1][k]+val[ii][j-1][k]+val[ii-1][j-1][k]));
            nzX = 0.5 * ( (val[ii][j][k+1]+val[ii][j][k+1]+val[ii-1][j][k+1])
                        - (val[ii][j][k-1]+val[ii][j][k-1]+val[ii-1][j][k-1]));
            normalize(nxX,nyX,nzX);
        
            nxY = 0.5 * ( (val[ii][j-1][k]+val[ii][j][k]+val[ii][j+1][k])
                        - (val[ii-1][j-1][k]+val[ii-1][j][k]+val[ii-1][j+1][k]));
            nyY = copysign(1.0,+(val[ii][j+1][k]-val[ii][j-1][k]));
            nzY = 0.5 * ( (val[ii][j-1][k+1]+val[ii][j][k+1]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j][k-1]+val[ii][j+1][k-1]));
            normalize(nxY,nyY,nzY);

            nxZ = 0.5 * ( (val[ii][j][k-1]+val[ii][j][k]+val[ii][j][k+1])
                        - (val[ii-1][j][k-1]+val[ii-1][j][k]+val[ii-1][j][k+1]));
            nyZ = 0.5 * ( (val[ii][j+1][k-1]+val[ii][j+1][k]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j-1][k]+val[ii][j-1][k+1]));
            nzZ = copysign(1.0,+(val[ii][j][k+1]-val[ii][j][k-1]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,ii,j,k);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;  // val[jj-1]=val[jj]
            nxX = copysign(1.0,+(val[i+1][jj][k]-val[i-1][jj][k]));
            nyX = 0.5 * ( (val[i+1][jj+1][k]+val[i][jj+1][k]+val[i-1][jj+1][k])
                        - (val[i+1][jj  ][k]+val[i][jj  ][k]+val[i-1][jj  ][k]));
            nzX = 0.5 * ( (val[i+1][jj][k+1]+val[i][jj][k+1]+val[i-1][jj][k+1])
                        - (val[i+1][jj][k-1]+val[i][jj][k-1]+val[i-1][jj][k-1]));
            normalize(nxX,nyX,nzX);

            nxY = 0.5 * ( (val[i+1][jj  ][k]+val[i+1][jj][k]+val[i+1][jj+1][k])
                        - (val[i-1][jj  ][k]+val[i-1][jj][k]+val[i-1][jj+1][k]));
            nyY = copysign(1.0,+(val[i][jj+1][k]-val[i][jj  ][k]));
            nzY = 0.5 * ( (val[i][jj  ][k+1]+val[i][jj][k+1]+val[i][jj+1][k+1])
                        - (val[i][jj  ][k-1]+val[i][jj][k-1]+val[i][jj+1][k-1]));
            normalize(nxY,nyY,nzY);

            nxZ = 0.5 * ( (val[i+1][jj][k-1]+val[i+1][jj][k]+val[i+1][jj][k+1])
                        - (val[i-1][jj][k-1]+val[i-1][jj][k]+val[i-1][jj][k+1]));
            nyZ = 0.5 * ( (val[i][jj+1][k-1]+val[i][jj+1][k]+val[i][jj+1][k+1])
                        - (val[i][jj  ][k-1]+val[i][jj  ][k]+val[i][jj  ][k+1]));
            nzZ = copysign(1.0,+(val[i][jj][k+1]-val[i][jj][k-1]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,jj,k);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1; // val[jj+1] = val[jj]
            nxX = copysign(1.0,+(val[i+1][jj][k]-val[i-1][jj][k]));
            nyX = 0.5 * ( (val[i+1][jj  ][k]+val[i][jj  ][k]+val[i-1][jj  ][k])
                        - (val[i+1][jj-1][k]+val[i][jj-1][k]+val[i-1][jj-1][k]));
            nzX = 0.5 * ( (val[i+1][jj][k+1]+val[i][jj][k+1]+val[i-1][jj][k+1])
                        - (val[i+1][jj][k-1]+val[i][jj][k-1]+val[i-1][jj][k-1]));
            normalize(nxX,nyX,nzX);

            nxY = 0.5 * ( (val[i+1][jj-1][k]+val[i+1][jj][k]+val[i+1][jj  ][k])
                        - (val[i-1][jj-1][k]+val[i-1][jj][k]+val[i-1][jj  ][k]));
            nyY = copysign(1.0,+(val[i][jj  ][k]-val[i][jj-1][k]));
            nzY = 0.5 * ( (val[i][jj-1][k+1]+val[i][jj][k+1]+val[i][jj  ][k+1])
                        - (val[i][jj-1][k-1]+val[i][jj][k-1]+val[i][jj  ][k-1]));
            normalize(nxY,nyY,nzY);

            nxZ = 0.5 * ( (val[i+1][jj][k-1]+val[i+1][jj][k]+val[i+1][jj][k+1])
                        - (val[i-1][jj][k-1]+val[i-1][jj][k]+val[i-1][jj][k+1]));
            nyZ = 0.5 * ( (val[i][jj  ][k-1]+val[i][jj  ][k]+val[i][jj  ][k+1])
                        - (val[i][jj-1][k-1]+val[i][jj-1][k]+val[i][jj-1][k+1]));
            nzZ = copysign(1.0,+(val[i][jj][k+1]-val[i][jj][k-1]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,jj,k);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1; //val[kk-1] = val[kk]
            nxX = copysign(1.0,+(val[i+1][j][kk]-val[i-1][j][kk]));
            nyX = 0.5 * ( (val[i+1][j+1][kk]+val[i][j+1][kk]+val[i-1][j+1][kk])
                        - (val[i+1][j-1][kk]+val[i][j-1][kk]+val[i-1][j-1][kk]));
            nzX = 0.5 * ( (val[i+1][j][kk+1]+val[i][j][kk+1]+val[i-1][j][kk+1])
                        - (val[i+1][j][kk  ]+val[i][j][kk  ]+val[i-1][j][kk  ]));
            normalize(nxX,nyX,nzX);

            nxY = 0.5 * ( (val[i+1][j-1][kk]+val[i+1][j][kk]+val[i+1][j+1][kk])
                        - (val[i-1][j-1][kk]+val[i-1][j][kk]+val[i-1][j+1][kk]));
            nyY = copysign(1.0,+(val[i][j+1][kk]-val[i][j-1][kk]));
            nzY = 0.5 * ( (val[i][j-1][kk+1]+val[i][j][kk+1]+val[i][j+1][kk+1])
                        - (val[i][j-1][kk  ]+val[i][j][kk  ]+val[i][j+1][kk  ]));
            normalize(nxY,nyY,nzY);

            nxZ = 0.5 * ( (val[i+1][j][kk  ]+val[i+1][j][kk]+val[i+1][j][kk+1])
                        - (val[i-1][j][kk  ]+val[i-1][j][kk]+val[i-1][j][kk+1]));
            nyZ = 0.5 * ( (val[i][j+1][kk  ]+val[i][j+1][kk]+val[i][j+1][kk+1])
                        - (val[i][j-1][kk  ]+val[i][j-1][kk]+val[i][j-1][kk+1]));
            nzZ = copysign(1.0,+(val[i][j][kk+1]-val[i][j][kk  ]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,kk);
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1; // val[kk+1] = val[kk]
            nxX = copysign(1.0,+(val[i+1][j][kk]-val[i-1][j][kk]));
            nyX = 0.5 * ( (val[i+1][j+1][kk]+val[i][j+1][kk]+val[i-1][j+1][kk])
                        - (val[i+1][j-1][kk]+val[i][j-1][kk]+val[i-1][j-1][kk]));
            nzX = 0.5 * ( (val[i+1][j][kk  ]+val[i][j][kk  ]+val[i-1][j][kk  ])
                        - (val[i+1][j][kk-1]+val[i][j][kk-1]+val[i-1][j][kk-1]));
            normalize(nxX,nyX,nzX);

            nxY = 0.5 * ( (val[i+1][j-1][kk]+val[i+1][j][kk]+val[i+1][j+1][kk])
                        - (val[i-1][j-1][kk]+val[i-1][j][kk]+val[i-1][j+1][kk]));
            nyY = copysign(1.0,+(val[i][j+1][kk]-val[i][j-1][kk]));
            nzY = 0.5 * ( (val[i][j-1][kk  ]+val[i][j][kk  ]+val[i][j+1][kk  ])
                        - (val[i][j-1][kk-1]+val[i][j][kk-1]+val[i][j+1][kk-1]));
            normalize(nxY,nyY,nzY);

            nxZ = 0.5 * ( (val[i+1][j][kk-1]+val[i+1][j][kk]+val[i+1][j][kk  ])
                        - (val[i-1][j][kk-1]+val[i-1][j][kk]+val[i-1][j][kk  ]));
            nyZ = 0.5 * ( (val[i][j+1][kk-1]+val[i][j+1][kk]+val[i][j+1][kk  ])
                        - (val[i][j-1][kk-1]+val[i][j-1][kk]+val[i][j-1][kk  ]));
            nzZ = copysign(1.0,+(val[i][j][kk  ]-val[i][j][kk-1]));
            normalize(nxZ,nyZ,nzZ);

            selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,kk);
          }
        }
      }
    }
#endif
  }

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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k]+val[i  ][j][k+1]));
      nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i  ][j][k-1]+val[i  ][j][k]));
      nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(nxY,nyY,nzY);

      nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k]+val[i][j  ][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
                  -(val[i][j-1][k  ]+val[i][j][k  ]));
      normalize(nxY,nyY,nzY);

      nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
                  -(val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
      normalize(nxX,nyX,nzX);

      nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(nxY,nyY,nzY);

      nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(nxZ,nyZ,nzZ);

      selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
               - (val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
               - (val[i][j][k-1]+val[i][j+1][k-1]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    nzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                -(val[i][j][k-1]+val[i][j+1][k-1]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    nyZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
    normalize(nxX,nyX,nzX);

    nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));
    normalize(nxY,nyY,nzY);

    nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
    normalize(nxZ,nyZ,nzZ);

    selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
  }
}


/******************************************************************************/
void VOF::selectMax(const real nxX, const real nyX, const real nzX,
               const real nxY, const real nyY, const real nzY,
               const real nxZ, const real nyZ, const real nzZ,
               const int i,    const int j,    const int k    ) {

  if (fabs(nxX)<fabs(nyY)) {
    if (fabs(nyY)<fabs(nzZ)) {
      nx[i][j][k]=nxZ;
      ny[i][j][k]=nyZ;
      nz[i][j][k]=nzZ;
    } else {
      nx[i][j][k]=nxY;
      ny[i][j][k]=nyY;
      nz[i][j][k]=nzY;
    }
  } else {
    if (fabs(nxX)<fabs(nzZ)) {
      nx[i][j][k]=nxZ;
      ny[i][j][k]=nyZ;
      nz[i][j][k]=nzZ;
    } else {
      nx[i][j][k]=nxX;
      ny[i][j][k]=nyX;
      nz[i][j][k]=nzX;
    }
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
  normalize(nxX,nyX,nzX);

  nxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(nxY,nyY,nzY);

  nxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
  normalize(nxX,nyX,nzX);

  nxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
             - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(nxY,nyY,nzY);

  nxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
  normalize(nxX,nyX,nzX);

  nxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
  nzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(nxY,nyY,nzY);

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
  normalize(nxX,nyX,nzX);

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]));
  nyY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
  nzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
              -(val[i][j-1][k-1]+val[i][j][k-1]));
  normalize(nxY,nyY,nzY);

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
  normalize(nxX,nyX,nzX);

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
  normalize(nxY,nyY,nzY);

  nxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k]+val[i-1][j][k+1]));
  nyZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k]+val[i][j-1][k+1]));
  nzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
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
  normalize(nxX,nyX,nzX);

  nxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  nyY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  nzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
              -(val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(nxY,nyY,nzY);

  nxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
              -(val[i-1][j][k-1]+val[i-1][j][k]));
  nyZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
              -(val[i][j-1][k-1]+val[i][j-1][k]));
  nzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
  normalize(nxZ,nyZ,nzZ);

  selectMax(nxX,nyX,nzX,nxY,nyY,nzY,nxZ,nyZ,nzZ,i,j,k);
}
