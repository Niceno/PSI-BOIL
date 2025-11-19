#include "model.h"

#include "../Plot/plot.h"

/******************************************************************************/
void Model::tau_wall(Momentum * mom, const Scalar & dist,
                     Vector * fbnd, Scalar * y_pl) const {

  static Scalar tau_x( * mom->domain() );
  static Scalar tau_y( * mom->domain() );
  static Scalar tau_z( * mom->domain() );
  tau_x = dist.shape();
  tau_y = dist.shape();
  tau_z = dist.shape();

  if( y_pl ) {
    assert(y_pl->domain() == mom->domain() );
    assert(y_pl->domain() == dist.domain() );
    * y_pl = 0.0;
  }

  /* take handy aliases */
  const Vector & uvw   =   mom->val();
  const Matter & fluid = * mom->fluid();
  const Body   & ibody =   mom->domain()->ibody(); 

  real tau_w,  tau_w_mean = 0.0;
  real y_plus, y_plus_mean = 0.0;
  int  cnt = 0;

  real tx, ty, tz;
  Comp m;

  /*--------+
  |         |
  |  y min  |
  |         |
  +--------*/
  if( dist.bc().type_here( Dir::jmin(), BndType::dirichlet()) ) {

    tau_x = 0.0;
    tau_y = 0.0;
    tau_z = 0.0;

    /* compute tau in scalar cells */
    const int j=dist.sj();
    for_vik(dist,i,k) 
      if( ibody.on(i,j,k) ) {
        wall_function( uvw, fluid, dist, & tau_w, & y_plus, i, j, k, &tx, &ty, &tz );
        y_plus_mean += y_plus; tau_w_mean += tau_w; cnt++;
        if(y_pl) (*y_pl)[i][j][k] = y_plus;

        tau_x[i][j][k] = dist.dSy(i,j,k) * tau_w * tx;
        tau_y[i][j][k] = dist.dSy(i,j,k) * tau_w * ty;
        tau_z[i][j][k] = dist.dSy(i,j,k) * tau_w * tz;
      } 
  }

  if( dist.bc().type( Dir::jmin(), BndType::dirichlet()) ) {
    tau_x.exchange();
    tau_y.exchange();
    tau_z.exchange();
  }

  if( dist.bc().type_here( Dir::jmin(), BndType::dirichlet()) ) {
    const int j=dist.sj();
    /* update momentum sources */
    for_vmik(uvw,Comp::u(),i,k) 
      if( tau_x[i][j][k] != 0.0 && tau_x[i-1][j][k] != 0.0 ) {
        mom->A[0]->c[i][j][k] -= mom->A[0]->s[i][j][k];
        mom->A[0]->s[i][j][k]  = 0.0;
        (*fbnd)[Comp::u()][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
      }
////for_vmik(uvw,Comp::v(),i,k)                                     // 
////  if( tau_y[i][j][k] != 0.0 ) {                                 //
////    mom->A[1]->c[i][j+1][k] -= mom->A[1]->s[i][j+1][k];         //  
////    mom->A[1]->s[i][j+1][k]  = 0.0;                             //
////    (*fbnd)[Comp::v()][i][j+1][k] -= 1.0 * (tau_y[i][j+1][k]);  //
////  }                                                             //
    for_vmik(uvw,Comp::w(),i,k) 
      if( tau_z[i][j][k] != 0.0 && tau_z[i][j][k-1] != 0.0 ) {
        mom->A[2]->c[i][j][k] -= mom->A[2]->s[i][j][k];
        mom->A[2]->s[i][j][k]  = 0.0;
        (*fbnd)[Comp::w()][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
      }
  }

  /*--------+
  |         |
  |  y max  |
  |         |
  +--------*/
  if( dist.bc().type_here( Dir::jmax(), BndType::dirichlet()) ) {

    tau_x = 0.0;
    tau_y = 0.0;
    tau_z = 0.0;

    /* compute tau in scalar cells */
    const int j=dist.ej();
    for_vik(dist,i,k) 
      if( ibody.on(i,j,k) ) {
        wall_function( uvw, fluid, dist, & tau_w, & y_plus, i, j, k, &tx, &ty, &tz );
        y_plus_mean += y_plus; tau_w_mean += tau_w; cnt++;
        if(y_pl) (*y_pl)[i][j][k] = y_plus;

        tau_x[i][j][k] = dist.dSy(i,j,k) * tau_w * tx;
        tau_y[i][j][k] = dist.dSy(i,j,k) * tau_w * ty;
        tau_z[i][j][k] = dist.dSy(i,j,k) * tau_w * tz;
      } 
  }

  if( dist.bc().type( Dir::jmax(), BndType::dirichlet()) ) {
    tau_x.exchange();
    tau_y.exchange();
    tau_z.exchange();
  }

  if( dist.bc().type_here( Dir::jmax(), BndType::dirichlet()) ) {
    const int j=dist.ej();
    /* update momentum sources */
    for_vmik(uvw,Comp::u(),i,k) 
      if( tau_x[i][j][k] != 0.0 && tau_x[i-1][j][k] != 0.0 ) {
        mom->A[0]->c[i][j][k] -= mom->A[0]->n[i][j][k];
        mom->A[0]->n[i][j][k]  = 0.0;
        (*fbnd)[Comp::u()][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
      }
////for_vmik(uvw,Comp::v(),i,k)                                     // 
////  if( tau_y[i][j][k] != 0.0 ) {                                 //
////    mom->A[1]->c[i][j][k] -= mom->A[1]->n[i][j][k];             // 
////    mom->A[1]->n[i][j][k]  = 0.0;                               //
////    (*fbnd)[Comp::v()][i][j][k] -= 1.0 * (tau_y[i][j][k]);      //
////  }                                                             //
    for_vmik(uvw,Comp::w(),i,k) 
      if( tau_z[i][j][k] != 0.0 && tau_z[i][j][k-1] != 0.0 ) {
        mom->A[2]->c[i][j][k] -= mom->A[2]->n[i][j][k];
        mom->A[2]->n[i][j][k]  = 0.0;
        (*fbnd)[Comp::w()][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
      }
  }

  /*--------+
  |         |
  |  k min  |
  |         |
  +--------*/
  if( dist.bc().type_here( Dir::kmin(), BndType::dirichlet()) ) {

    tau_x = 0.0;
    tau_y = 0.0;
    tau_z = 0.0;

    /* compute tau in scalar cells */
    const int k=dist.sk();
    for_vij(dist,i,j) 
      if( ibody.on(i,j,k) ) {
        wall_function( uvw, fluid, dist, & tau_w, & y_plus, i, j, k, &tx, &ty, &tz );
        y_plus_mean += y_plus; tau_w_mean += tau_w; cnt++;
        if(y_pl) (*y_pl)[i][j][k] = y_plus;

        tau_x[i][j][k] = dist.dSz(i,j,k) * tau_w * tx;
        tau_y[i][j][k] = dist.dSz(i,j,k) * tau_w * ty;
        tau_z[i][j][k] = dist.dSz(i,j,k) * tau_w * tz;
      } 
  }

  if( dist.bc().type( Dir::kmin(), BndType::dirichlet()) ) {
    tau_x.exchange();
    tau_y.exchange();
    tau_z.exchange();
  }
    
  if( dist.bc().type_here( Dir::kmin(), BndType::dirichlet()) ) {
    const int k=dist.sk();
    /* update momentum sources */
    for_vmij(uvw,Comp::u(),i,j) 
      if( tau_x[i][j][k] != 0.0 && tau_x[i-1][j][k] != 0.0 ) {
        mom->A[0]->c[i][j][k] -= mom->A[0]->b[i][j][k];
        mom->A[0]->b[i][j][k]  = 0.0;
        (*fbnd)[Comp::u()][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
      }
    for_vmij(uvw,Comp::v(),i,j) 
      if( tau_y[i][j][k] != 0.0 && tau_y[i][j-1][k] != 0.0 ) {
        mom->A[1]->c[i][j][k] -= mom->A[1]->b[i][j][k];
        mom->A[1]->b[i][j][k]  = 0.0;
        (*fbnd)[Comp::v()][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
      }
////for_vmij(uvw,Comp::w(),i,j)                                     //
////  if( tau_z[i][j][k] != 0.0 ) {                                 //
////    mom->A[2]->c[i][j][k+1] -= mom->A[2]->b[i][j][k+1];         // 
////    mom->A[2]->b[i][j][k+1]  = 0.0;                             //
////    (*fbnd)[Comp::w()][i][j][k+1] -= 0.5 * (tau_z[i][j][k]);    //
////  }                                                             //
  }

  /*--------+
  |         |
  |  k max  |
  |         |
  +--------*/
  if( dist.bc().type_here( Dir::kmax(), BndType::dirichlet()) ) {

    tau_x = 0.0;
    tau_y = 0.0;
    tau_z = 0.0;

    /* compute tau in scalar cells */
    const int k=dist.ek();
    for_vij(dist,i,j) 
      if( ibody.on(i,j,k) ) {
        wall_function( uvw, fluid, dist, & tau_w, & y_plus, i, j, k, &tx, &ty, &tz );
        y_plus_mean += y_plus; tau_w_mean += tau_w; cnt++;
        if(y_pl) (*y_pl)[i][j][k] = y_plus;

        tau_x[i][j][k] = dist.dSz(i,j,k) * tau_w * tx;
        tau_y[i][j][k] = dist.dSz(i,j,k) * tau_w * ty;
        tau_z[i][j][k] = dist.dSz(i,j,k) * tau_w * tz;
      } 
  }

  if( dist.bc().type( Dir::kmax(), BndType::dirichlet()) ) {
    tau_x.exchange();
    tau_y.exchange();
    tau_z.exchange();
  }
    
  if( dist.bc().type_here( Dir::kmax(), BndType::dirichlet()) ) {
    const int k=dist.ek();
    /* update momentum sources */
    for_vmij(uvw,Comp::u(),i,j) 
      if( tau_x[i][j][k] != 0.0 && tau_x[i-1][j][k] != 0.0 ) {
        mom->A[0]->c[i][j][k] -= mom->A[0]->t[i][j][k];
        mom->A[0]->t[i][j][k]  = 0.0;
        (*fbnd)[Comp::u()][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
      }
    for_vmij(uvw,Comp::v(),i,j) 
      if( tau_y[i][j][k] != 0.0 && tau_y[i][j-1][k] != 0.0 ) {
        mom->A[1]->c[i][j][k] -= mom->A[1]->t[i][j][k];
        mom->A[1]->t[i][j][k]  = 0.0;
        (*fbnd)[Comp::v()][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
      }
////for_vmij(uvw,Comp::w(),i,j)                                     //
////  if( tau_z[i][j][k] != 0.0 ) {                                 // 
////    mom->A[2]->c[i][j][k] -= mom->A[2]->t[i][j][k];             // 
////    mom->A[2]->t[i][j][k]  = 0.0;                               //    
////    (*fbnd)[Comp::w()][i][j][k] -= (tau_z[i][j][k]);            //
////  }                                                             //
  }

  /*----------------+
  |                 |
  |  immersed body  |
  |                 |
  +----------------*/
  /* create offset indices */
  int io[6]={-1,+1, 0, 0, 0, 0};
  int jo[6]={ 0, 0,-1,+1, 0, 0};
  int ko[6]={ 0, 0, 0, 0,-1,+1};

  if( ibody.ncall() > 0 ) {

    tau_x = 0.0;
    tau_y = 0.0;
    tau_z = 0.0;
  //if( ibody.nccells() > 0 ) {
    //OMS(gotcha!);
    for_vijk(dist,i,j,k) {
      if( ibody.on(i,j,k) ) {
        for(int dir=0; dir<6; dir++) {
          if( ibody.off(i+io[dir], j+jo[dir], k+ko[dir]) ) {

            wall_function( uvw, fluid, dist, & tau_w, & y_plus
                         , i, j, k, &tx, &ty, &tz );
            y_plus_mean += y_plus; 
            tau_w_mean  += tau_w; 
            cnt++;
            if(y_pl) (*y_pl)[i][j][k] = y_plus;

            real dS;
            if( dir == 0 || dir == 1 ) dS = dist.dSx(i,j,k);
            if( dir == 2 || dir == 3 ) dS = dist.dSy(i,j,k);
            if( dir == 4 || dir == 5 ) dS = dist.dSz(i,j,k);

            tau_x[i][j][k] += dS * tau_w * tx;
            tau_y[i][j][k] += dS * tau_w * ty;
            tau_z[i][j][k] += dS * tau_w * tz;
          }
        }
      }
    } 

    tau_x.exchange();
    tau_y.exchange();
    tau_z.exchange();
    
    /*--------------------------+
    |  update momentum sources  | -> not quite finished yet
    +--------------------------*/

    m = Comp::u();
    for_vmijk(uvw,m,i,j,k) {
      if( ibody.on(m,i,j,k) ) {
        if( ibody.off(m,i,j-1,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->s[i][j][k];
          mom->A[~m]->s[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
        } 
        if( ibody.off(m,i,j+1,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->n[i][j][k];
          mom->A[~m]->n[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
        } 
        if( ibody.off(m,i,j,k-1) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->b[i][j][k];
          mom->A[~m]->b[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
        } 
        if( ibody.off(m,i,j,k+1) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->t[i][j][k];
          mom->A[~m]->t[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_x[i][j][k] + tau_x[i-1][j][k]);
        } 
      }
    }

    m = Comp::v();
    for_vmijk(uvw,m,i,j,k) {
      if( ibody.on(m,i,j,k) ) {
        if( ibody.off(m,i-1,j,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->w[i][j][k];
          mom->A[~m]->w[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
        } 
        if( ibody.off(m,i+1,j,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->e[i][j][k];
          mom->A[~m]->e[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
        } 
        if( ibody.off(m,i,j,k-1) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->b[i][j][k];
          mom->A[~m]->b[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
        } 
        if( ibody.off(m,i,j,k+1) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->t[i][j][k];
          mom->A[~m]->t[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_y[i][j][k] + tau_y[i][j-1][k]);
        } 
      }
    }

    m = Comp::w();
    for_vmijk(uvw,m,i,j,k) {
      if( ibody.on(m,i,j,k) ) {
        if( ibody.off(m,i-1,j,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->w[i][j][k];
          mom->A[~m]->w[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
        } 
        if( ibody.off(m,i+1,j,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->e[i][j][k];
          mom->A[~m]->e[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
        } 
        if( ibody.off(m,i,j-1,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->s[i][j][k];
          mom->A[~m]->s[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
        } 
        if( ibody.off(m,i,j+1,k) ) {
          mom->A[~m]->c[i][j][k] -= mom->A[~m]->n[i][j][k];
          mom->A[~m]->n[i][j][k]  = 0.0;
          (*fbnd)[m][i][j][k] -= 0.5 * (tau_z[i][j][k] + tau_z[i][j][k-1]);
        } 
      }
    }
  }

  boil::cart.sum_real(&y_plus_mean);
  boil::cart.sum_real(&tau_w_mean );
  boil::cart.sum_int (&cnt        );

  if( cnt > 0 ) {
    y_plus_mean /= (real)cnt;
    tau_w_mean  /= (real)cnt;
    boil::oout<<"model_tau_wall:cnt,y_plus,tau_w= "<<cnt<<" "<<y_plus_mean<<" "
              <<tau_w_mean<<"\n";
  }
  else
    OMS(No near wall cells found);

  if(y_pl) y_pl->exchange();

  //debug: boil::plot->plot(tau_x, tau_y, tau_z, "tau");
}
