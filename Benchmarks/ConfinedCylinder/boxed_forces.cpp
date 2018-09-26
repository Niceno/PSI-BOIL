/******************************************************************************/
real viscous_force_x(const Vector & uvw, int I, int J, int K,
                     const Range<int> & ir,  
                     const Range<int> & jr) {

  real vf = 0.0;

  for(int i=ir.first(); i<=ir.last(); i++) {
    /* south */
     {
      const int j = jr.first();
      const real u_1  = uvw[Comp::u()][I+i][J+j]  [K];
      const real u_1s = uvw[Comp::u()][I+i][J+j-1][K];
      vf += MU * (u_1s - u_1) * uvw.dSy(Comp::u(),I+i,J+j,K) 
                              / uvw.dys(Comp::u(),J+j);    /* [N] */
     }
    /* north */
     {
      const int j = jr.last();
      const real u_1  = uvw[Comp::u()][I+i][J+j]  [K];
      const real u_1n = uvw[Comp::u()][I+i][J+j+1][K];
      vf += MU * (u_1n - u_1) * uvw.dSy(Comp::u(),I+i,J+j,K) 
                              / uvw.dyn(Comp::u(),J+j);    /* [N] */
     }
  }

  for(int j=jr.first(); j<=jr.last(); j++) {
    /* west */
     {
      const int i = ir.first();
      const real u_1  = uvw[Comp::u()][I+i]  [J+j][K];
      const real u_1w = uvw[Comp::u()][I+i-1][J+j][K];
      vf += MU * (u_1w - u_1) * uvw.dSx(Comp::u(),I+i,J+j,K) 
                              / uvw.dxw(Comp::u(),I+i);    /* [N] */ 
     }
    /* east */
     {
      const int i = ir.last();
      const real u_1  = uvw[Comp::u()][I+i]  [J+j][K];
      const real u_1e = uvw[Comp::u()][I+i+1][J+j][K];
      vf += MU * (u_1e - u_1) * uvw.dSx(Comp::u(),I+i,J+j,K) 
                              / uvw.dxe(Comp::u(),I+i);    /* [N] */
     }
  }

  /* return normalized value */
  return vf / uvw.dzc(Comp::u(), K);
}

/******************************************************************************/
real viscous_force_y(const Vector & uvw, int I, int J, int K,
                     const Range<int> & ir,  
                     const Range<int> & jr) {

  real vf = 0.0;

  for(int i=ir.first(); i<=ir.last(); i++) {
    /* south */
     {
      const int j = jr.first();
      const real v_1  = uvw[Comp::v()][I+i][J+j]  [K];
      const real v_1s = uvw[Comp::v()][I+i][J+j-1][K];
      vf += MU * (v_1s - v_1) * uvw.dSy(Comp::v(),I+i,J+j,K) 
                              / uvw.dys(Comp::v(),J+j);    /* [N] */
     }
    /* north */
     {
      const int j = jr.last();
      const real v_1  = uvw[Comp::v()][I+i][J+j]  [K];
      const real v_1n = uvw[Comp::v()][I+i][J+j+1][K];
      vf += MU * (v_1n - v_1) * uvw.dSy(Comp::v(),I+i,J+j,K) 
                              / uvw.dyn(Comp::v(),J+j);    /* [N] */
     }
  }

  for(int j=jr.first(); j<=jr.last(); j++) {
    /* west */
     {
      const int i = ir.first();
      const real v_1  = uvw[Comp::v()][I+i]  [J+j][K];
      const real v_1w = uvw[Comp::v()][I+i-1][J+j][K];
      vf += MU * (v_1w - v_1) * uvw.dSx(Comp::v(),I+i,J+j,K) 
                              / uvw.dxw(Comp::v(),I+i);    /* [N] */ 
     }
    /* east */
     {
      const int i = ir.last();
      const real v_1  = uvw[Comp::v()][I+i]  [J+j][K];
      const real v_1e = uvw[Comp::v()][I+i+1][J+j][K];
      vf += MU * (v_1e - v_1) * uvw.dSx(Comp::v(),I+i,J+j,K) 
                              / uvw.dxe(Comp::v(),I+i);    /* [N] */
     }
  }

  /* return normalized value */
  return vf / uvw.dzc(Comp::v(), K);
}

/******************************************************************************/
real convective_force_x(const Vector & uvw, int I, int J, int K,   
                        const Range<int> & ir,  
                        const Range<int> & jr) {

  real cf = 0.0;

  for(int i=ir.first(); i<=ir.last(); i++) {
    /* south */
     {
      const int j = jr.first();
      const real u_1s = 0.5 * ( uvw[Comp::u()][I+i]  [J+j-1][K] + 
                                uvw[Comp::u()][I+i]  [J+j]  [K]);
      const real v_1s = 0.5 * ( uvw[Comp::v()][I+i]  [J+j]  [K] +
                                uvw[Comp::v()][I+i-1][J+j]  [K] );
      cf += RHO * u_1s * v_1s * uvw.dSy(Comp::u(),I+i,J+j,K);
     }
    /* north */
     {
      const int j = jr.last();
      const real u_1n = 0.5 * ( uvw[Comp::u()][I+i]  [J+j+1][K] + 
                                uvw[Comp::u()][I+i]  [J+j]  [K]);
      const real v_1n = 0.5 * ( uvw[Comp::v()][I+i]  [J+j+1][K] +
                                uvw[Comp::v()][I+i-1][J+j+1][K] );
      cf -= RHO * u_1n * v_1n * uvw.dSy(Comp::u(),I+i,J+j,K);
     }
  }

  for(int j=jr.first(); j<=jr.last(); j++) {
    /* west */
     {
      const int i = ir.first();
      const real u_1  = uvw[Comp::u()][I+i]  [J+j][K];
      const real u_1w = uvw[Comp::u()][I+i-1][J+j][K];
      cf += RHO * 0.25*(u_1w+u_1)*(u_1w+u_1)*uvw.dSx(Comp::u(),I+i,J+j,K);
     }
    /* east */
     {
      const int i = ir.last();
      const real u_1  = uvw[Comp::u()][I+i]  [J+j][K];
      const real u_1e = uvw[Comp::u()][I+i+1][J+j][K];
      cf -= RHO * 0.25*(u_1e+u_1)*(u_1e+u_1)*uvw.dSx(Comp::u(),I+i,J+j,K);
     }
  }

  /* return normalized value */
  return cf / uvw.dzc(Comp::u(), K);
}

/******************************************************************************/
real convective_force_y(const Vector & uvw, int I, int J, int K,   
                        const Range<int> & ir,  
                        const Range<int> & jr) {

  real cf = 0.0;

  for(int i=ir.first(); i<=ir.last(); i++) {
    /* south */
     {
      const int j = jr.first();
      const real v_1  = uvw[Comp::v()][I+i][J+j]  [K];
      const real v_1s = uvw[Comp::v()][I+i][J+j-1][K];
      cf += RHO * 0.25*(v_1s+v_1)*(v_1s+v_1)*uvw.dSy(Comp::v(),I+i,J+j,K);
     }
    /* north */
     {
      const int j = jr.last();
      const real v_1  = uvw[Comp::v()][I+i][J+j]  [K];
      const real v_1n = uvw[Comp::v()][I+i][J+j+1][K];
      cf -= RHO * 0.25*(v_1n+v_1)*(v_1n+v_1)*uvw.dSy(Comp::v(),I+i,J+j,K);
     }
  }

  for(int j=jr.first(); j<=jr.last(); j++) {
    /* west */
     {
      const int i = ir.first();
      const real u_1w = 0.5 * ( uvw[Comp::u()][I+i]  [J+j-1][K] + 
                                uvw[Comp::u()][I+i]  [J+j]  [K]);
      const real v_1w = 0.5 * ( uvw[Comp::v()][I+i]  [J+j]  [K] +
                                uvw[Comp::v()][I+i-1][J+j]  [K] );
      cf += RHO * u_1w * v_1w * uvw.dSx(Comp::v(),I+i,J+j,K);
     }
    /* east */
     {
      const int i = ir.last();
      const real u_1e = 0.5 * ( uvw[Comp::u()][I+i+1][J+j-1][K] + 
                                uvw[Comp::u()][I+i+1][J+j]  [K]);
      const real v_1e = 0.5 * ( uvw[Comp::v()][I+i]  [J+j]  [K] +
                                uvw[Comp::v()][I+i+1][J+j]  [K] );
      cf -= RHO * u_1e * v_1e * uvw.dSx(Comp::v(),I+i,J+j,K);
     }
  }

  /* return normalized value */
  return cf / uvw.dzc(Comp::v(), K);
}

/******************************************************************************/
real innertia_x(const Vector & uvw, int I, int J, int K,
                const Range<int> & ir,
                const Range<int> & jr) {

  real mf = 0.0;

  const Body & b = uvw.domain()->ibody();

  for(int i=ir.first(); i<=ir.last(); i++) 
    for(int j=jr.first(); j<=jr.last(); j++) {
      if( b.fV(Comp::u(),I+i,J+j,K) > 0.5 ) {
        const real u_1 = uvw[Comp::u()][I+i][J+j][K];
        const real fV  = b.fV(Comp::u(),I+i,J+j,K);
        mf += u_1 * RHO * uvw.dV(Comp::u(),I+i,J+j,K) * fV; /* [kg m / s] */
      }
    }

  /* return normalized value */
  return mf / uvw.dzc(Comp::u(), K);
}

/******************************************************************************/
real innertia_y(const Vector & uvw, int I, int J, int K,
                const Range<int> & ir,
                const Range<int> & jr) {

  real mf = 0.0;

  const Body & b = uvw.domain()->ibody();

  for(int i=ir.first(); i<=ir.last(); i++) 
    for(int j=jr.first(); j<=jr.last(); j++) {
      if( b.fV(Comp::v(),I+i,J+j,K) > 0.5 ) {
        const real v_1 = uvw[Comp::v()][I+i][J+j][K];
        const real fV  = b.fV(Comp::v(),I+i,J+j,K);
        mf += v_1 * RHO * uvw.dV(Comp::v(),I+i,J+j,K) * fV; /* [kg m / s] */
      }
    }

  /* return normalized value */
  return mf / uvw.dzc(Comp::v(), K);
}

/******************************************************************************/
real pressure_force_x(const Scalar & p, int I, int J, int K,
                      const Range<int> & ir,  
                      const Range<int> & jr) {

  real pf = 0.0;

  for(int j=jr.first(); j<=jr.last(); j++) {
    /* west */
     {
      const int i = ir.first()-1;
      pf += p[I+i][J+j][K] * p.dSx(I+i,J+j,K); 
     }
    /* east */
     {
      const int i = ir.last();
      pf -= p[I+i][J+j][K] * p.dSx(I+i,J+j,K); 
     }
  }

  /* return normalized value */
  return pf / p.dzc(K);
}

/******************************************************************************/
real pressure_force_y(const Scalar & p, int I, int J, int K,
                      const Range<int> & ir,  
                      const Range<int> & jr) {

  real pf = 0.0;

  for(int i=ir.first(); i<=ir.last(); i++) {
    /* south */
     {
      const int j = jr.first()-1;
      pf += p[I+i][J+j][K] * p.dSy(I+i,J+j,K); 
     }
    /* north */
     {
      const int j = jr.last();
      pf -= p[I+i][J+j][K] * p.dSy(I+i,J+j,K); 
     }
  }

  /* return normalized value */
  return pf / p.dzc(K);
}

