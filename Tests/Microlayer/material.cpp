/******************************************************************************/
/* ------------ material properties */
  const real Mv = IF97::get_MW();
  const real muv = IF97::viscvap_p(prs);
  const real rhov = IF97::rhovap_p(prs);
  const real cpv = IF97::cpvap_p(prs)*rhov;
  const real lambdav = IF97::tcondvap_p(prs);

  const real mul = IF97::viscliq_p(prs);
  const real rhol = IF97::rholiq_p(prs);
  const real cpl = IF97::cpliq_p(prs)*rhol;
  const real lambdal = IF97::tcondliq_p(prs);

  const real sig = IF97::sigma97(tsat0_K);
  const real latent=IF97::hvap_p(prs)-IF97::hliq_p(prs);

  //const real betal = (7.03+(tsat0_K-273.15-100.)/160.*(22.1-7.03))*1e-4; /* roughly */

  /* estimate betal */
  const real betal = -(IF97::rholiq_p(1.01*prs)-IF97::rholiq_p(0.99*prs))
                     /(IF97::Tsat97(1.01*prs)-IF97::Tsat97(0.99*prs))/rhol;

  const real betav = 1./tsat0_K; /* ideal gas approximation */

  const real alpl = lambdal/cpl;
  const real alpv = lambdav/cpv;

  boil::oout<<"properties at pressure "<<prs<<boil::endl;
  boil::oout<<"vapprop= "<<Mv<<" "<<muv<<" "<<rhov<<" "
                         <<cpv/rhov<<" "<<lambdav<<" "<<betav<<" "<<alpv<<boil::endl;
  boil::oout<<"liqprop= "<<Mv<<" "<<mul<<" "<<rhol<<" "
                         <<cpl/rhol<<" "<<lambdal<<" "<<betal<<" "<<alpl<<boil::endl;
  boil::oout<<"twoprop= "<<tsat0_K<<" "<<sig<<" "<<latent<<boil::endl;

  const real Jal = cpl*deltat_nucl/(latent*rhov);
  boil::oout << "Jal= "<<Jal<<boil::endl;

  /* heater */
  /* sapphire */
  const real rhosol = 3980.0;
  const real trefsol = tsat0_K;
  //real cpsol = 0.9161e3; //D. A. Ditmars, et. al., J. Res. Nat. Bur. Stand., 87, (2), 159-163 (1982).
  real cpsol = 0.929e3;
  cpsol *= rhosol;
  const real lambdasol = 25.12;
  const real alpsol = lambdasol/cpsol;
#if 0
  /* ito */
  const real rhoheat = 7140;
  const real cpheat = 2.58e6;
  const real lambdaheat = 10.2;
#else
  /* titanium */
  const real rhoheat = 4506.;
  const real cpheat = 544.3*rhoheat;
  const real lambdaheat = 17.;
#endif

  /* ghost distance */
  /* HT resistance */
  real resistance_liq(0.0);
  real accomult = 2.*accommodation/(2.-accommodation);
  if(accommodation>0.0) {
    resistance_liq = Schrage::calculate_heat_transfer_resistance(tsat0_K,rhov,Mv,latent)
                   /(accomult*0.5);
  } 

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  TwoLevelMatter vapor(d), liquid(d);
  vapor  .mu    (muv);
  vapor  .rho   (rhov);
  vapor  .cp    (cpv);  /* J/m3 */
  vapor  .lambda(lambdav);
  vapor  .mmass (Mv);
  vapor  .beta  (betav);
  liquid.mu    (mul);
  liquid.rho   (rhol);
  liquid.cp    (cpl);   /* J/m3 */
  liquid.lambda(lambdal);
  liquid.mmass (Mv);
  liquid.beta  (betal);

  TwoLevelMatter mixed(liquid, vapor, & c);
  mixed.sigma(sig);
  mixed.latent(latent);

  TwoLevelMatter substrate(d), heater(d);
  substrate.rho    (rhosol);
  substrate.cp     (cpsol);
  substrate.lambda (lambdasol);
  heater.rho    (rhoheat);
  heater.cp     (cpheat);
  heater.lambda (lambdaheat);
  
  TwoLevelMatter solid(substrate,heater,&csub);
