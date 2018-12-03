void compare_areas(Scalar & phi, Heaviside & pc, Scalar & heaviadens) {
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

}
