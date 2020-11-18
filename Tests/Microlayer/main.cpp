#include "Include/psi-boil.h"
#include <fstream>
#include <iostream>
#include <fenv.h>
#include <iterator>
#define _GNU_SOURCE 1
#if 1
static void __attribute__ ((constructor)) trapfpe(void)
{
  /* Enable some exceptions. At startup all exceptions are masked. */
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}
#endif

#include "header.cpp"
#include "aux.cpp"

/******************************************************************************/
int main(int argc, char ** argv) {

  boil::timer.start();

  if(argc<3){
    boil::oout<<"Two command line arguments required!"<<"\n";
    boil::oout<<"./Boil wmin file-with-inputs"<<"\n";
    boil::oout<<"Inputs:\n";
    boil::oout<<"- wall superheat [K]\n";
    boil::oout<<"- outlet superheat [K]\n";
    boil::oout<<"- nucleation superheat [K]\n";
    boil::oout<<"- heat flux [W/m2], will be converted to vol. source\n";
    boil::oout<<"- cangle [deg]\n";
    boil::oout<<"- accommodation coefficient [-] \n";
    boil::oout<<"- pressure [atm]\n";
    boil::oout<<"- initradius [um]\n";
    boil::oout<<"- LZsolid [um]\n";
    boil::oout<<"- LZheater[um], 0. indicates a Dirac heater\n";
    boil::oout<<"- multiplicative length in sol [-]\n";
    boil::oout<<"- multiplicative length in x [-]\n";
    boil::oout<<"- multiplicative length in z [-]\n";
    boil::oout<<"- max cell aspect ratio [-]\n";
    boil::oout<<"- NZsolid [.int.]\n";
    boil::oout<<"- NXtot [.int.]\n";
    boil::oout<<"- NZtot [.int.]\n";
    boil::oout<<"- transitional no. cells in sol [.int.]\n";
    boil::oout<<"- transitional no. cells in X [.int.]\n";
    boil::oout<<"- transitional no. cells in Z [.int.]\n";
    boil::oout<<"- case flag [.int.], indicates initialisation and time loop\n";

    exit(0);
  }

/* case flags: 
 * 0 - single-phase simulation for temperature distribution
 * 1 - multi-phase, on coarse mesh, uniform temperatures
 * 2 - multi-phase, on coarse mesh, temperatures from a file
 * 3 - multi-phase, on fine mesh, uniform temperatures
 * 4 - multi-phase, on fine mesh, temperatures from a file
 */

/******************************************************************************/
/* ------------ input from command line */
  int wmin=atoi(argv[1]);
  boil::oout<<"wmin= "<<wmin<<"\n";

  const std::string deck(argv[2]);
  boil::oout<<"Reading input file << "<<deck<<" >>\n";

  /* important: set all values to zero to start with */
  real surftens_dt_coef(0.);
  real deltat_wall(0.), deltat_out(0.), deltat_nucl(0.), qflux(0.);
  real cangle(0.), accommodation(0.), prs(0.), radius(0.);
  real LZsol(0.), LZheat(0.);
  real LXmult(0.), LZmult(0.), LZsolmult(0.);
  real AR(0.);
  int NZsol(0), NXtot(0), NZtot(0);
  int NZsol_trans(0), NX_trans(0), NZ_trans(0);
  int factor_x(0), factor_z(0);
  int case_flag(0);

  std::vector<real*> readreal({&surftens_dt_coef,
                               &deltat_wall,
                               &deltat_out,
                               &deltat_nucl,
                               &qflux,
                               &cangle,
                               &accommodation,
                               &prs,
                               &radius,
                               &LZsol,
                               &LZheat,
                               &LZsolmult,
                               &LXmult,
                               &LZmult,
                               &AR
                              });

  std::vector<int*> readint({&NZsol,
                             &NXtot,
                             &NZtot,
                             &NZsol_trans,
                             &NX_trans,
                             &NZ_trans,
                             &factor_x,
                             &factor_z,
                             &case_flag
                            });

  /* only root reads */
  int allread(0);
  std::ifstream input;
  if(boil::cart.iam()==0) {
    input.open(deck,std::ios::in);
    if(!input.fail()) {
      for(auto & rr : readreal)
        allread += read_variable(input,*rr);
      for(auto & ri : readint)
        allread += read_variable(input,*ri);
    } else {
      boil::oout<<"Input file not found. Exiting."<<"\n";
      exit(0);
    }
    input.close();
  }

  /* post-process */
  prs *= 101325.; /* to pascal */
  radius *= 1e-6; /* to m */
  LZsol *= 1e-6; /* to m */
  LZheat *= 1e-6; /* to m */

  boil::cart.sum_int(&allread);
  if(allread) {
    boil::oout<<"Not all variables read. Exiting."<<"\n";
    exit(0);
  } else {
    boil::oout<<"All variables read successfully."<<"\n";
  }

  /* distribute read values over all processors (bit of an mpi trick) */
  for(auto & rr : readreal) {
    boil::cart.sum_real(rr);
  }
  for(auto & ri : readint) {
    boil::cart.sum_int(ri);
  }

#include "settings.cpp"

#include "domain.cpp"

#include "variables.cpp"

#include "material.cpp"

#include "solvers.cpp"

  /*-------------------+
  |  check if restart  |
  +-------------------*/
  boil::test_irun();
  boil::set_irun(1);

#include "init.cpp"

#ifndef SETUP_ONLY
  if(case_flag==0) {
  #include "timeloop_singlephase.cpp"
  } else if(case_flag<3) {
  #include "timeloop_multiphase_coarse.cpp"
  } else {
    OMS(Underdevelopment. Exiting.);
    exit(0);
  }
#endif

  boil::oout << "Finished." << boil::endl;
  boil::timer.stop();
  boil::timer.report();

  return 0;
}

