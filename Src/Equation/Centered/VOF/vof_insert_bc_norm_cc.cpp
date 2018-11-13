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

        real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;

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

    if( val.bc().type(b) == BndType::symmetry() ) {

      /*-----------+
      |  Symmetry  |
      +-----------*/

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      if(d != Dir::undefined()) {

        real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;

        if(d == Dir::imin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i+1; // val[ii-1]=val[ii]
            mxX = copysign(1.0,+(val[ii+1][j][k]-val[ii][j][k]));
            myX = 0.5 * ( (val[ii+1][j+1][k]+val[ii][j+1][k]+val[ii][j+1][k])
                        - (val[ii+1][j-1][k]+val[ii][j-1][k]+val[ii][j-1][k]));
            mzX = 0.5 * ( (val[ii+1][j][k+1]+val[ii][j][k+1]+val[ii][j][k+1])
                        - (val[ii+1][j][k-1]+val[ii][j][k-1]+val[ii][j][k-1]));
            normalize(mxX,myX,mzX);

            mxY = 1.0 * ( (val[ii+1][j-1][k]+val[ii+1][j][k]+val[ii+1][j+1][k])
                        - (val[ii  ][j-1][k]+val[ii  ][j][k]+val[ii  ][j+1][k]));
            myY = copysign(1.0,+(val[ii][j+1][k]-val[ii][j-1][k]));
            mzY = 0.5 * ( (val[ii][j-1][k+1]+val[ii][j][k+1]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j][k-1]+val[ii][j+1][k-1]));
            normalize(mxY,myY,mzY);

            mxZ = 1.0 * ( (val[ii+1][j][k-1]+val[ii+1][j][k]+val[ii+1][j][k+1])
                        - (val[ii  ][j][k-1]+val[ii  ][j][k]+val[ii  ][j][k+1]));
            myZ = 0.5 * ( (val[ii][j+1][k-1]+val[ii][j+1][k]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j-1][k]+val[ii][j-1][k+1]));
            mzZ = copysign(1.0,+(val[ii][j][k+1]-val[ii][j][k-1]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,ii,j,k);
          }
        }
        if(d == Dir::imax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int ii=i-1;  // val[ii+1]=val[ii]
            mxX = copysign(1.0,+(val[ii][j][k]-val[ii-1][j][k]));
            myX = 0.5 * ( (val[ii][j+1][k]+val[ii][j+1][k]+val[ii-1][j+1][k])
                        - (val[ii][j-1][k]+val[ii][j-1][k]+val[ii-1][j-1][k]));
            mzX = 0.5 * ( (val[ii][j][k+1]+val[ii][j][k+1]+val[ii-1][j][k+1])
                        - (val[ii][j][k-1]+val[ii][j][k-1]+val[ii-1][j][k-1]));
            normalize(mxX,myX,mzX);
        
            mxY = 0.5 * ( (val[ii][j-1][k]+val[ii][j][k]+val[ii][j+1][k])
                        - (val[ii-1][j-1][k]+val[ii-1][j][k]+val[ii-1][j+1][k]));
            myY = copysign(1.0,+(val[ii][j+1][k]-val[ii][j-1][k]));
            mzY = 0.5 * ( (val[ii][j-1][k+1]+val[ii][j][k+1]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j][k-1]+val[ii][j+1][k-1]));
            normalize(mxY,myY,mzY);

            mxZ = 0.5 * ( (val[ii][j][k-1]+val[ii][j][k]+val[ii][j][k+1])
                        - (val[ii-1][j][k-1]+val[ii-1][j][k]+val[ii-1][j][k+1]));
            myZ = 0.5 * ( (val[ii][j+1][k-1]+val[ii][j+1][k]+val[ii][j+1][k+1])
                        - (val[ii][j-1][k-1]+val[ii][j-1][k]+val[ii][j-1][k+1]));
            mzZ = copysign(1.0,+(val[ii][j][k+1]-val[ii][j][k-1]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,ii,j,k);
          }
        }
        if(d == Dir::jmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j+1;  // val[jj-1]=val[jj]
            mxX = copysign(1.0,+(val[i+1][jj][k]-val[i-1][jj][k]));
            myX = 0.5 * ( (val[i+1][jj+1][k]+val[i][jj+1][k]+val[i-1][jj+1][k])
                        - (val[i+1][jj  ][k]+val[i][jj  ][k]+val[i-1][jj  ][k]));
            mzX = 0.5 * ( (val[i+1][jj][k+1]+val[i][jj][k+1]+val[i-1][jj][k+1])
                        - (val[i+1][jj][k-1]+val[i][jj][k-1]+val[i-1][jj][k-1]));
            normalize(mxX,myX,mzX);

            mxY = 0.5 * ( (val[i+1][jj  ][k]+val[i+1][jj][k]+val[i+1][jj+1][k])
                        - (val[i-1][jj  ][k]+val[i-1][jj][k]+val[i-1][jj+1][k]));
            myY = copysign(1.0,+(val[i][jj+1][k]-val[i][jj  ][k]));
            mzY = 0.5 * ( (val[i][jj  ][k+1]+val[i][jj][k+1]+val[i][jj+1][k+1])
                        - (val[i][jj  ][k-1]+val[i][jj][k-1]+val[i][jj+1][k-1]));
            normalize(mxY,myY,mzY);

            mxZ = 0.5 * ( (val[i+1][jj][k-1]+val[i+1][jj][k]+val[i+1][jj][k+1])
                        - (val[i-1][jj][k-1]+val[i-1][jj][k]+val[i-1][jj][k+1]));
            myZ = 0.5 * ( (val[i][jj+1][k-1]+val[i][jj+1][k]+val[i][jj+1][k+1])
                        - (val[i][jj  ][k-1]+val[i][jj  ][k]+val[i][jj  ][k+1]));
            mzZ = copysign(1.0,+(val[i][jj][k+1]-val[i][jj][k-1]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,jj,k);
          }
        }
        if(d == Dir::jmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int jj=j-1; // val[jj+1] = val[jj]
            mxX = copysign(1.0,+(val[i+1][jj][k]-val[i-1][jj][k]));
            myX = 0.5 * ( (val[i+1][jj  ][k]+val[i][jj  ][k]+val[i-1][jj  ][k])
                        - (val[i+1][jj-1][k]+val[i][jj-1][k]+val[i-1][jj-1][k]));
            mzX = 0.5 * ( (val[i+1][jj][k+1]+val[i][jj][k+1]+val[i-1][jj][k+1])
                        - (val[i+1][jj][k-1]+val[i][jj][k-1]+val[i-1][jj][k-1]));
            normalize(mxX,myX,mzX);

            mxY = 0.5 * ( (val[i+1][jj-1][k]+val[i+1][jj][k]+val[i+1][jj  ][k])
                        - (val[i-1][jj-1][k]+val[i-1][jj][k]+val[i-1][jj  ][k]));
            myY = copysign(1.0,+(val[i][jj  ][k]-val[i][jj-1][k]));
            mzY = 0.5 * ( (val[i][jj-1][k+1]+val[i][jj][k+1]+val[i][jj  ][k+1])
                        - (val[i][jj-1][k-1]+val[i][jj][k-1]+val[i][jj  ][k-1]));
            normalize(mxY,myY,mzY);

            mxZ = 0.5 * ( (val[i+1][jj][k-1]+val[i+1][jj][k]+val[i+1][jj][k+1])
                        - (val[i-1][jj][k-1]+val[i-1][jj][k]+val[i-1][jj][k+1]));
            myZ = 0.5 * ( (val[i][jj  ][k-1]+val[i][jj  ][k]+val[i][jj  ][k+1])
                        - (val[i][jj-1][k-1]+val[i][jj-1][k]+val[i][jj-1][k+1]));
            mzZ = copysign(1.0,+(val[i][jj][k+1]-val[i][jj][k-1]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,jj,k);
          }
        }
        if(d == Dir::kmin()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k+1; //val[kk-1] = val[kk]
            mxX = copysign(1.0,+(val[i+1][j][kk]-val[i-1][j][kk]));
            myX = 0.5 * ( (val[i+1][j+1][kk]+val[i][j+1][kk]+val[i-1][j+1][kk])
                        - (val[i+1][j-1][kk]+val[i][j-1][kk]+val[i-1][j-1][kk]));
            mzX = 0.5 * ( (val[i+1][j][kk+1]+val[i][j][kk+1]+val[i-1][j][kk+1])
                        - (val[i+1][j][kk  ]+val[i][j][kk  ]+val[i-1][j][kk  ]));
            normalize(mxX,myX,mzX);

            mxY = 0.5 * ( (val[i+1][j-1][kk]+val[i+1][j][kk]+val[i+1][j+1][kk])
                        - (val[i-1][j-1][kk]+val[i-1][j][kk]+val[i-1][j+1][kk]));
            myY = copysign(1.0,+(val[i][j+1][kk]-val[i][j-1][kk]));
            mzY = 0.5 * ( (val[i][j-1][kk+1]+val[i][j][kk+1]+val[i][j+1][kk+1])
                        - (val[i][j-1][kk  ]+val[i][j][kk  ]+val[i][j+1][kk  ]));
            normalize(mxY,myY,mzY);

            mxZ = 0.5 * ( (val[i+1][j][kk  ]+val[i+1][j][kk]+val[i+1][j][kk+1])
                        - (val[i-1][j][kk  ]+val[i-1][j][kk]+val[i-1][j][kk+1]));
            myZ = 0.5 * ( (val[i][j+1][kk  ]+val[i][j+1][kk]+val[i][j+1][kk+1])
                        - (val[i][j-1][kk  ]+val[i][j-1][kk]+val[i][j-1][kk+1]));
            mzZ = copysign(1.0,+(val[i][j][kk+1]-val[i][j][kk  ]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,kk);
          }
        }
        if(d == Dir::kmax()){
          for_vijk( val.bc().at(b), i,j,k ){
            int kk=k-1; // val[kk+1] = val[kk]
            mxX = copysign(1.0,+(val[i+1][j][kk]-val[i-1][j][kk]));
            myX = 0.5 * ( (val[i+1][j+1][kk]+val[i][j+1][kk]+val[i-1][j+1][kk])
                        - (val[i+1][j-1][kk]+val[i][j-1][kk]+val[i-1][j-1][kk]));
            mzX = 0.5 * ( (val[i+1][j][kk  ]+val[i][j][kk  ]+val[i-1][j][kk  ])
                        - (val[i+1][j][kk-1]+val[i][j][kk-1]+val[i-1][j][kk-1]));
            normalize(mxX,myX,mzX);

            mxY = 0.5 * ( (val[i+1][j-1][kk]+val[i+1][j][kk]+val[i+1][j+1][kk])
                        - (val[i-1][j-1][kk]+val[i-1][j][kk]+val[i-1][j+1][kk]));
            myY = copysign(1.0,+(val[i][j+1][kk]-val[i][j-1][kk]));
            mzY = 0.5 * ( (val[i][j-1][kk  ]+val[i][j][kk  ]+val[i][j+1][kk  ])
                        - (val[i][j-1][kk-1]+val[i][j][kk-1]+val[i][j+1][kk-1]));
            normalize(mxY,myY,mzY);

            mxZ = 0.5 * ( (val[i+1][j][kk-1]+val[i+1][j][kk]+val[i+1][j][kk  ])
                        - (val[i-1][j][kk-1]+val[i-1][j][kk]+val[i-1][j][kk  ]));
            myZ = 0.5 * ( (val[i][j+1][kk-1]+val[i][j+1][kk]+val[i][j+1][kk  ])
                        - (val[i][j-1][kk-1]+val[i][j-1][kk]+val[i][j-1][kk  ]));
            mzZ = copysign(1.0,+(val[i][j][kk  ]-val[i][j][kk-1]));
            normalize(mxZ,myZ,mzZ);

            selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,kk);
          }
        }
      }
    }
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
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]));
      mzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j][k]+val[i  ][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      mzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-min & j-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    for_k(k) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      mzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]));
      myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      mzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
      myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-min & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int k=sk();
    for_j(j) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      myX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i  ][j][k]+val[i  ][j][k+1]));
      myZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-min & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int k=ek();
    for_j(j) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
      myX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]));
      mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i  ][j][k-1]+val[i  ][j][k]));
      myZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-max & j-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    for_k(k) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j  ][k]+val[i-1][j  ][k]));
      mzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      mzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-max & j-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    for_k(k) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      mzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]));
      myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      mzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
                 - (val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
      myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-max & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int k=sk();
    for_j(j) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      myX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      mzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i][j][k  ]+val[i-1][j][k  ]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
                 - (val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      myZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line i-max & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int k=ek();
    for_j(j) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      myX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j-1][k]+val[i-1][j-1][k]));
      mzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
                 - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
                 - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      myZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line j-min & k-min */
  if(val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=sj();
    int k=sk();
    for_i(i) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
      mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));
      normalize(mxX,myX,mzX);

      mxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      mzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                  -(val[i][j][k  ]+val[i][j+1][k  ]));
      normalize(mxY,myY,mzY);

      mxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      myZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                  -(val[i][j  ][k]+val[i][j  ][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line j-min & k-max */
  if(val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=sj();
    int k=ek();
    for_i(i) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
      mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      mzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line j-max & k-min */
  if(val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int j=ej();
    int k=sk();
    for_i(i) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
      mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
                  -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));
      normalize(mxX,myX,mzX);

      mxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
                  -(val[i][j-1][k  ]+val[i][j][k  ]));
      normalize(mxY,myY,mzY);

      mxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
                  -(val[i-1][j][k]+val[i-1][j][k+1]));
      myZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                  -(val[i][j-1][k]+val[i][j-1][k+1]));
      mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* line j-max & k-max */
  if(val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int j=ej();
    int k=ek();
    for_i(i) {
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
                  -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
      mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
                  -(val[i-1][j-1][k]+val[i-1][j][k]));
      myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
      mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
                  -(val[i][j-1][k-1]+val[i][j][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                  -(val[i][j-1][k-1]+val[i][j-1][k]));
      mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
    }
  }

  /* corner i-max & j-min & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=sk();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    myX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    mzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    mzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    myZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-max & j-min & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=ek();
      real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
      mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
      myX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                  -(val[i][j  ][k]+val[i-1][j  ][k]));
      mzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                  -(val[i][j][k-1]+val[i-1][j][k-1]));
      normalize(mxX,myX,mzX);

      mxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                  -(val[i-1][j][k]+val[i-1][j+1][k]));
      myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
      mzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                  -(val[i][j][k-1]+val[i][j+1][k-1]));
      normalize(mxY,myY,mzY);

      mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                  -(val[i-1][j][k-1]+val[i-1][j][k]));
      myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                  -(val[i][j  ][k-1]+val[i][j  ][k]));
      mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
      normalize(mxZ,myZ,mzZ);

      selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-min & j-min & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int j=sj();
    int k=sk();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                -(val[i+1][j  ][k]+val[i][j  ][k]));
    mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                -(val[i+1][j][k  ]+val[i][j][k  ]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    mzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
               - (val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    myZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-min & j-min & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int j=sj();
    int k=ek();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k])
                -(val[i+1][j  ][k]+val[i][j  ][k]));
    mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                -(val[i+1][j][k-1]+val[i][j][k-1]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i+1][j][k]+val[i+1][j+1][k])
                -(val[i  ][j][k]+val[i  ][j+1][k]));
    myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    mzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
               - (val[i][j][k-1]+val[i][j+1][k-1]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-min & j-max & k-min */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    int k=sk();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                -(val[i+1][j-1][k]+val[i][j-1][k]));
    mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1])
                -(val[i+1][j][k  ]+val[i][j][k  ]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i+1][j][k]+val[i+1][j][k+1])
                -(val[i  ][j][k]+val[i  ][j][k+1]));
    myZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-min & j-max & k-max */
  if(val.bc().type_here(Dir::imin(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=si();
    int j=ej();
    int k=ek();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
    myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k])
                -(val[i+1][j-1][k]+val[i][j-1][k]));
    mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ])
                -(val[i+1][j][k-1]+val[i][j][k-1]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k])
                -(val[i  ][j-1][k]+val[i  ][j][k]));
    myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k])
                -(val[i  ][j][k-1]+val[i  ][j][k]));
    myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-max & j-min & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=sk();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    myX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    mzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    mzY = 1.0 * ((val[i][j][k+1]+val[i][j+1][k+1])
                -(val[i][j][k  ]+val[i][j+1][k  ]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    myZ = 1.0 * ((val[i][j+1][k]+val[i][j+1][k+1])
                -(val[i][j  ][k]+val[i][j  ][k+1]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-max & j-min & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmin(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=sj();
    int k=ek();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    myX = 1.0 * ((val[i][j+1][k]+val[i-1][j+1][k])
                -(val[i][j  ][k]+val[i-1][j  ][k]));
    mzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                -(val[i][j][k-1]+val[i-1][j][k-1]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i  ][j][k]+val[i  ][j+1][k])
                -(val[i-1][j][k]+val[i-1][j+1][k]));
    myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
    mzY = 1.0 * ((val[i][j][k  ]+val[i][j+1][k  ])
                -(val[i][j][k-1]+val[i][j+1][k-1]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k])
                -(val[i][j  ][k-1]+val[i][j  ][k]));
    mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-max & j-max & k-min */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmin(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    int k=sk();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    myX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                -(val[i][j-1][k]+val[i-1][j-1][k]));
    mzX = 1.0 * ((val[i][j][k+1]+val[i-1][j][k+1])
                -(val[i][j][k  ]+val[i-1][j][k  ]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1])
               - (val[i][j-1][k  ]+val[i][j][k  ]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i  ][j][k]+val[i  ][j][k+1])
                -(val[i-1][j][k]+val[i-1][j][k+1]));
    myZ = 1.0 * ((val[i][j  ][k]+val[i][j  ][k+1])
                -(val[i][j-1][k]+val[i][j-1][k+1]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }

  /* corner i-max & j-max & k-max */
  if(val.bc().type_here(Dir::imax(), BndType::wall()) &&
     val.bc().type_here(Dir::jmax(), BndType::wall()) &&
     val.bc().type_here(Dir::kmax(), BndType::wall())   ) {
    int i=ei();
    int j=ej();
    int k=ek();
    real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
    mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
    myX = 1.0 * ((val[i][j  ][k]+val[i-1][j  ][k])
                -(val[i][j-1][k]+val[i-1][j-1][k]));
    mzX = 1.0 * ((val[i][j][k  ]+val[i-1][j][k  ])
                -(val[i][j][k-1]+val[i-1][j][k-1]));
    normalize(mxX,myX,mzX);

    mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k])
               - (val[i-1][j-1][k]+val[i-1][j][k]));
    myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
    mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ])
               - (val[i][j-1][k-1]+val[i][j][k-1]));
    normalize(mxY,myY,mzY);

    mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k])
                -(val[i-1][j][k-1]+val[i-1][j][k]));
    myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k])
                -(val[i][j-1][k-1]+val[i][j-1][k]));
    mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
    normalize(mxZ,myZ,mzZ);

    selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
  }
}


/******************************************************************************/
void VOF::selectMax(const real mxX, const real myX, const real mzX,
               const real mxY, const real myY, const real mzY,
               const real mxZ, const real myZ, const real mzZ,
               const int i,    const int j,    const int k    ) {

  if (fabs(mxX)<fabs(myY)) {
    if (fabs(myY)<fabs(mzZ)) {
      nx[i][j][k]=mxZ;
      ny[i][j][k]=myZ;
      nz[i][j][k]=mzZ;
    } else {
      nx[i][j][k]=mxY;
      ny[i][j][k]=myY;
      nz[i][j][k]=mzY;
    }
  } else {
    if (fabs(mxX)<fabs(mzZ)) {
      nx[i][j][k]=mxZ;
      ny[i][j][k]=myZ;
      nz[i][j][k]=mzZ;
    } else {
      nx[i][j][k]=mxX;
      ny[i][j][k]=myX;
      nz[i][j][k]=mzX;
    }
  }
}

/******************************************************************************/
void VOF::norm_cc_imin(const Scalar & val,
                       const int i,const int j, const int k) {
  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i+1][j][k]-val[i][j][k]));
  myX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]));
  mzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]));
  normalize(mxX,myX,mzX);

  mxY = 1.0 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k]));
  myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  mzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(mxY,myY,mzY);

  mxZ = 1.0 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1]));
  myZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

/******************************************************************************/
void VOF::norm_cc_imax(const Scalar & val,
                       const int i,const int j, const int k) {
  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i][j][k]-val[i-1][j][k]));
  myX = 0.5 * ((val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i][j-1][k]+val[i-1][j-1][k]));
  mzX = 0.5 * ((val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i][j][k-1]+val[i-1][j][k-1]));
  normalize(mxX,myX,mzX);

  mxY = 1.0 * ((val[i  ][j-1][k]+val[i  ][j][k]+val[i  ][j+1][k])
             - (val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  mzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
             - (val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(mxY,myY,mzY);

  mxZ = 1.0 * ((val[i  ][j][k-1]+val[i  ][j][k]+val[i  ][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  myZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

/******************************************************************************/
void VOF::norm_cc_jmin(const Scalar & val,
                       const int i,const int j, const int k) {
  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  myX = 1.0 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k]));
  mzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));
  normalize(mxX,myX,mzX);

  mxY = 0.5 * ((val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j][k]+val[i-1][j+1][k]));
  myY = copysign(1.0,+(val[i][j+1][k]-val[i][j][k]));
  mzY = 0.5 * ((val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(mxY,myY,mzY);

  mxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  myZ = 1.0 * ((val[i][j+1][k-1]+val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1]));
  mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

/******************************************************************************/
void VOF::norm_cc_jmax(const Scalar & val,
                       const int i,const int j, const int k) {
  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  myX = 1.0 * ((val[i+1][j  ][k]+val[i][j  ][k]+val[i-1][j  ][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  mzX = 0.5 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));
  normalize(mxX,myX,mzX);

  mxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]));
  myY = copysign(1.0,+(val[i][j][k]-val[i][j-1][k]));
  mzY = 0.5 * ((val[i][j-1][k+1]+val[i][j][k+1])
              -(val[i][j-1][k-1]+val[i][j][k-1]));
  normalize(mxY,myY,mzY);

  mxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k-1]+val[i-1][j][k]+val[i-1][j][k+1]));
  myZ = 1.0 * ((val[i][j  ][k-1]+val[i][j  ][k]+val[i][j  ][k+1])
              -(val[i][j-1][k-1]+val[i][j-1][k]+val[i][j-1][k+1]));
  mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k-1]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

/******************************************************************************/
void VOF::norm_cc_kmin(const Scalar & val,
                       const int i,const int j, const int k) {
  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  myX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  mzX = 1.0 * ((val[i+1][j][k+1]+val[i][j][k+1]+val[i-1][j][k+1])
              -(val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ]));
  normalize(mxX,myX,mzX);

  mxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  mzY = 1.0 * ((val[i][j-1][k+1]+val[i][j][k+1]+val[i][j+1][k+1])
              -(val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ]));
  normalize(mxY,myY,mzY);

  mxZ = 0.5 * ((val[i+1][j][k]+val[i+1][j][k+1])
              -(val[i-1][j][k]+val[i-1][j][k+1]));
  myZ = 0.5 * ((val[i][j+1][k]+val[i][j+1][k+1])
              -(val[i][j-1][k]+val[i][j-1][k+1]));
  mzZ = copysign(1.0,+(val[i][j][k+1]-val[i][j][k]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

/******************************************************************************/
void VOF::norm_cc_kmax(const Scalar & val,
                       const int i,const int j, const int k) {

  real mxX, myX, mzX, mxY, myY, mzY, mxZ, myZ, mzZ;
  mxX = copysign(1.0,+(val[i+1][j][k]-val[i-1][j][k]));
  myX = 0.5 * ((val[i+1][j+1][k]+val[i][j+1][k]+val[i-1][j+1][k])
              -(val[i+1][j-1][k]+val[i][j-1][k]+val[i-1][j-1][k]));
  mzX = 1.0 * ((val[i+1][j][k  ]+val[i][j][k  ]+val[i-1][j][k  ])
              -(val[i+1][j][k-1]+val[i][j][k-1]+val[i-1][j][k-1]));
  normalize(mxX,myX,mzX);

  mxY = 0.5 * ((val[i+1][j-1][k]+val[i+1][j][k]+val[i+1][j+1][k])
              -(val[i-1][j-1][k]+val[i-1][j][k]+val[i-1][j+1][k]));
  myY = copysign(1.0,+(val[i][j+1][k]-val[i][j-1][k]));
  mzY = 1.0 * ((val[i][j-1][k  ]+val[i][j][k  ]+val[i][j+1][k  ])
              -(val[i][j-1][k-1]+val[i][j][k-1]+val[i][j+1][k-1]));
  normalize(mxY,myY,mzY);

  mxZ = 0.5 * ((val[i+1][j][k-1]+val[i+1][j][k])
              -(val[i-1][j][k-1]+val[i-1][j][k]));
  myZ = 0.5 * ((val[i][j+1][k-1]+val[i][j+1][k])
              -(val[i][j-1][k-1]+val[i][j-1][k]));
  mzZ = copysign(1.0,+(val[i][j][k]-val[i][j][k-1]));
  normalize(mxZ,myZ,mzZ);

  selectMax(mxX,myX,mzX,mxY,myY,mzY,mxZ,myZ,mzZ,i,j,k);
}

