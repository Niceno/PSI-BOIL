  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor_1(dom_1), liquid_1(dom_1), sapphire_1(dom_1);
  vapor_1  .mu    (1.255e-5);
  vapor_1  .rho   (0.597);
  vapor_1  .cp    (2030*0.597);
  vapor_1  .lambda(0.025);
  liquid_1.mu    (0.28e-3);
  liquid_1.rho   (958.4);
  liquid_1.cp    (4215.9*958.4);
  liquid_1.lambda(0.679);
  sapphire_1.rho    (3980.0);
  sapphire_1.cp     (750.0*3980.0);
  sapphire_1.lambda (35.0);
  Matter mixed_1(liquid_1, vapor_1, & step_1);
  mixed_1.sigma(5.9e-2);

