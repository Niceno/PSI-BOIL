void compare_areas(Scalar & phi, 
                   Heaviside & pc, Scalar & heaviadens,
                   Scalar & gradadens,
                   Scalar & scaledA,
                   Scalar & scaledB) {
 
  boil::timer.start("MC adens");

  /* marching cube */
  real MCarea(0.0);
  int MCcount(0);
  for_vijk(phi,i,j,k) {
    real mcareaplus = pc.area(i,j,k);
    heaviadens[i][j][k] = mcareaplus/heaviadens.dV(i,j,k);
    if (mcareaplus > boil::atto) MCcount++;
    MCarea += mcareaplus;
  }
  
  boil::oout<<"MC: "<<MCcount<<" "<<MCarea<<boil::endl;
   
  boil::timer.stop("MC adens");
  boil::timer.start("GC adens");

  /* gradclr */
  real GCarea(0.0);
  int GCcount(0);
  real SAarea(0.0);
  int SAcount(0);
  real SBarea(0.0);
  int SBcount(0);
  for_vijk(phi,i,j,k) {

    real gradient = 0.0;
    real dxm = -phi.dxw(i);
    real dxp = phi.dxe(i);
    real dym = -phi.dys(j);
    real dyp = phi.dyn(j);
    real dzm = -phi.dzb(k);
    real dzp = phi.dzt(k);

    real c00 = phi[i][j][k];
    real cxm = phi[i-1][j][k];
    real cxp = phi[i+1][j][k];
    real cym = phi[i][j-1][k];
    real cyp = phi[i][j+1][k];
    real czm = phi[i][j][k-1];
    real czp = phi[i][j][k+1];

    real gradx = (cxp-cxm)/(dxp-dxm); 
    real grady = (cyp-cym)/(dyp-dym); 
    real gradz = (czp-czm)/(dzp-dzm); 

    gradient = sqrt(gradx*gradx+grady*grady+gradz*gradz);

    /* non-scaled */
    gradadens[i][j][k] = gradient;
    if (gradadens[i][j][k] > boil::atto) GCcount++;
    GCarea += gradadens[i][j][k] * phi.dV(i,j,k);

    /* scaled A = 2*alpha*gradient */
    scaledA[i][j][k] = gradient*c00*2.0;
    if (scaledA[i][j][k] > boil::atto) SAcount++;
    SAarea += scaledA[i][j][k] * phi.dV(i,j,k);

    /* scaled B = 6*alpha*(1-alpha)*gradient */
    scaledB[i][j][k] = gradient*c00*(1.-c00)*6.0;
    if (scaledB[i][j][k] > boil::atto) SBcount++;
    SBarea += scaledB[i][j][k] * phi.dV(i,j,k);
  }

  boil::oout<<"GC: "<<GCcount<<" "<<GCarea<<boil::endl;
  boil::oout<<"SA: "<<SAcount<<" "<<SAarea<<boil::endl;
  boil::oout<<"SB: "<<SBcount<<" "<<SBarea<<boil::endl;

  boil::timer.stop("GC adens");

}
