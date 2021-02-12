/******************************************************************************/
/* ------------ material properties */
#if CASE == 0
  const real tsat0_K = IF97::Tsat97(prs);
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

  /* estimate betal */
  const real betal = -(IF97::rholiq_p(1.01*prs)-IF97::rholiq_p(0.99*prs))
                     /(IF97::Tsat97(1.01*prs)-IF97::Tsat97(0.99*prs))/rhol;

  boil::oout<<"properties at pressure "<<prs<<boil::endl;
#elif CASE == 1
  /* FC-72 from Cao(2019) and FC-72 3M product sheet */
  const real tsat0_K = 55.7+273.15;                      
  const real Mv = 338e-3;
  const real muv = 1.2e-5;
  const real rhov = 13.33;
  const real cpv = 894.*rhov;
  const real lambdav = 0.0129;

  const real mul = 4.4e-4;
  const real rhol = 1740.-2.61*(tsat0_K-273.15);
  const real cpl = (1014+1.554*(tsat0_K-273.15))*rhol;
  const real lambdal = 0.06-0.00011*(tsat0_K-273.15);

  const real sig = 7.9e-3;
  const real latent = 76900.;

  const real betal = 0.00156;

  /* fit-correlation from Raj(2012) */
  cangle = 9.37*std::pow(deltat_nucl,0.54);

  boil::oout<<"properties for FC-72 (atmospheric)"<<boil::endl;
#elif CASE == 2
  /* ethanol from Carey ISBN 978-1-59169-035-1 */
  const real tsat0_K = 78.37+273.15;
  const real Mv = 46.07e-3;
  const real muv = 10.4e-6;
  const real rhov = 1.435;
  const real cpv = 1.83e3*rhov;
  const real lambdav = 19.9e-3;

  const real mul = 428.7e-6;
  const real rhol = 757.0;
  const real cpl = 3.0e3*rhol;
  const real lambdal = 153.6e-3;

  const real sig = 1.77e-2;
  const real latent = 963.0e3;

  /* val at 50 degC, from Sun et al. (1988), doi: 10.1080/00319108808078584 */
  const real betal = 0.001151; 

  boil::oout<<"properties for ethanol (atmospheric)"<<boil::endl;
#elif CASE == 3
  /* sodium */
  boil::oout<<"Underdevelopment. Exiting."<<boil::endl;
  exit(0);
#elif CASE == 4
  /* prototype */
  const real tsat0_K = 400.;
  const real Mv = 40e-3;

  const real muv = 1e-5;
  const real rhov = 1.;
  const real cpv = 2e3*rhov;
  const real lambdav = 1e-2;

  const real mul = 1.25e-4;
  const real rhol = 1000.0;
  const real cpl = 4.0e3*rhol;
  const real lambdal = mul*cpl/rhol/freeinp;//5e-2;

  const real sig = 1e-2;
  const real latent = 1e6;

  const real betal = 1e-3;

  boil::oout<<"properties for prototype fluid"<<boil::endl;
#else
  boil::oout<<"Underdevelopment. Exiting."<<boil::endl;
#endif

  const real betav = 1./tsat0_K; /* ideal gas approximation */

  const real alpl = lambdal/cpl;
  const real alpv = lambdav/cpv;

  boil::oout<<"vapprop= "<<Mv<<" "<<muv<<" "<<rhov<<" "
                         <<cpv/rhov<<" "<<lambdav<<" "<<betav<<" "<<alpv<<boil::endl;
  boil::oout<<"liqprop= "<<Mv<<" "<<mul<<" "<<rhol<<" "
                         <<cpl/rhol<<" "<<lambdal<<" "<<betal<<" "<<alpl<<boil::endl;
  boil::oout<<"twoprop= "<<tsat0_K<<" "<<sig<<" "<<latent<<boil::endl;

  const real Jal = cpl*deltat_nucl/(latent*rhov);
  boil::oout << "Jal= "<<Jal<<" "<<Jal*sqrt(alpl)<<boil::endl;

  const real Prl = (mul/rhol)/alpl;
  boil::oout << "Prl= "<<Prl<<boil::endl;

  /* heater */
#ifndef SAKASHITA
  /* sapphire */
  const real rhosol = 3980.0;
  //real cpsol = 0.9161e3; //D. A. Ditmars, et. al., J. Res. Nat. Bur. Stand., 87, (2), 159-163 (1982).
  real cpsol = 0.929e3;
  const real lambdasol = 25.12;
#else
  /* nickel */ 
  const real rhosol = 8908.0;
  real cpsol = 444.0;
  const real lambdasol = 90.9;
#endif
  cpsol *= rhosol;
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

  /* HT resistance */
  real resistance_liq(0.0);
  if(accommodation>0.0) {
    resistance_liq = 
      TIF::calculate_heat_transfer_resistance
                     (tsat0_K,rhov,Mv,latent,accommodation);
  } 

  /* evaporative capillary number */
  real nanoscale = 10e-9;
  real Cae(0.0), thetahock(0.0);
  if(resistance_liq>0.) {
    Cae = mul*deltat_nucl/rhol/latent/resistance_liq/sig;
    thetahock = std::pow(12*Cae*log(0.5*dxmin/nanoscale),0.25);
  }
  boil::oout<<"Cae= "<<Cae<<" "<<thetahock<<boil::endl;

  /*----------------------+
  |  physical properties  |
  +----------------------*/
  Matter vapor(d), liquid(d);
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

  Matter mixed(liquid, vapor, & c);
  mixed.sigma(sig);
  mixed.latent(latent);

  Matter substrate(d), heater(d);
  substrate.rho    (rhosol);
  substrate.cp     (cpsol);
  substrate.lambda (lambdasol);
  heater.rho    (rhoheat);
  heater.cp     (cpheat);
  heater.lambda (lambdaheat);
  
  Matter solid(substrate,heater,&csub);
