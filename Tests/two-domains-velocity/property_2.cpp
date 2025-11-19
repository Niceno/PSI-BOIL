  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor_2(dom_2), liquid_2(dom_2), sapphire_2(dom_2);
  vapor_2  .mu    (1.255e-5);
  vapor_2  .rho   (0.597);
  vapor_2  .cp    (2030*0.597);
  vapor_2  .lambda(0.025);
  liquid_2.mu    (0.28e-3);
  liquid_2.rho   (958.4);
  liquid_2.cp    (4215.9*958.4);
  liquid_2.lambda(0.679);
  sapphire_2.rho    (3980.0);
  sapphire_2.cp     (750.0*3980.0);
  sapphire_2.lambda (35.0);
  Matter mixed_2(liquid_2, vapor_2, & step_2);
  mixed_2.sigma(5.9e-2);

