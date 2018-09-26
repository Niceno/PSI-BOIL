#ifndef VARIABLES_H
#define VARIABLES_H

#include <vector>

#include "structures.h"

/* global variable definitions */
extern std::vector<struct ele> elem;
extern std::vector<struct sid> side;
extern std::vector<struct nod> node;
extern std::vector<struct bc>  bound;

extern real dt, visc, CFL;  /* what is this for ? */ 
extern char name[256]; 
extern int  len;
extern std::vector<real> SxV, SyV, SxD, SyD;

extern class matrix Auvw, Aphi;
extern std::vector<real> u, v, uo, vo, 
                         bu_b,   bv_b,   /* boundary conditions */
                         bu_Co,  bv_Co,  /* convective terms */
                         bu_Coo, bv_Coo, /* convective terms */
                         bu_Do,  bv_Do,  /* diffusuve terms */
                         bu_n,   bv_n,   /* nonstationary terms */
                         bu_P,   bv_P,   /* pressure terms */
                         bu,     bv,     /* total fluxes */
                         phi, bphi, p, Flux, 
                         dp_dx, dp_dy;

#endif

/*-----------------------------------------------------------------------------+
 '$Id: variables.h,v 1.2 2015/09/01 09:51:02 niceno Exp $'/
+-----------------------------------------------------------------------------*/
