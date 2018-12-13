#include "../../Ravioli/comp.h"

#ifndef BROWSING_VECTOR_H
#define BROWSING_VECTOR_H

//--------------------
// through components 
//--------------------
#define for_m(m) for(Comp m=Comp::u(); m<=Comp::w(); m++)

// CHECK THESE BROWSERS IN PARALLEL!!!

//----
// 1D
//----
#define for_mi(m,i) for(int i=si(m); i<=ei(m); i++)
#define for_mj(m,j) for(int j=sj(m); j<=ej(m); j++)
#define for_mk(m,k) for(int k=sk(m); k<=ek(m); k++)

/* a */
#define for_ami(m,i) for(int i=0; i<ni(m); i++)
#define for_amj(m,j) for(int j=0; j<nj(m); j++)
#define for_amk(m,k) for(int k=0; k<nk(m); k++)

/* v */
#define for_vmi(v,m,i) for(int i=v.si(m); i<=v.ei(m); i++)
#define for_vmj(v,m,j) for(int j=v.sj(m); j<=v.ej(m); j++)
#define for_vmk(v,m,k) for(int k=v.sk(m); k<=v.ek(m); k++)

/* av */
#define for_avmi(v,m,i) for(int i=0; i<v.ni(m); i++)
#define for_avmj(v,m,j) for(int j=0; j<v.nj(m); j++)
#define for_avmk(v,m,k) for(int k=0; k<v.nk(m); k++)

#if 1

/* w = v which does not skip walls */
#define for_wmi(m,i) for(int i=1; i<ni(m)-1; i++)
#define for_wmj(m,j) for(int j=1; j<nj(m)-1; j++)
#define for_wmk(m,k) for(int k=1; k<nk(m)-1; k++)

/* wv */
#define for_wvmi(v,m,i) for(int i=1; i<v.ni(m)-1; i++)
#define for_wvmj(v,m,j) for(int j=1; j<v.nj(m)-1; j++)
#define for_wvmk(v,m,k) for(int k=1; k<v.nk(m)-1; k++)

#endif

//----
// 2D
//----
#define for_mij(m,i,j) for_mi(m,i) for_mj(m,j)
#define for_mik(m,i,k) for_mi(m,i) for_mk(m,k)
#define for_mjk(m,j,k) for_mj(m,j) for_mk(m,k)

/* a */
#define for_amij(m,i,j) for_ami(m,i) for_amj(m,j)
#define for_amik(m,i,k) for_ami(m,i) for_amk(m,k)
#define for_amjk(m,j,k) for_amj(m,j) for_amk(m,k)

/* v */
#define for_vmij(v,m,i,j) for_vmi(v,m,i) for_vmj(v,m,j)
#define for_vmik(v,m,i,k) for_vmi(v,m,i) for_vmk(v,m,k)
#define for_vmjk(v,m,j,k) for_vmj(v,m,j) for_vmk(v,m,k)

/* av */
#define for_avmij(v,m,i,j) for_avmi(v,m,i) for_avmj(v,m,j)
#define for_avmik(v,m,i,k) for_avmi(v,m,i) for_avmk(v,m,k)
#define for_avmjk(v,m,j,k) for_avmj(v,m,j) for_avmk(v,m,k)

#if 1
/* w */
#define for_wmij(m,i,j) for_wmi(m,i) for_wmj(m,j)
#define for_wmik(m,i,k) for_wmi(m,i) for_wmk(m,k)
#define for_wmjk(m,j,k) for_wmj(m,j) for_wmk(m,k)

/* wv */
#define for_wvmij(v,m,i,j) for_wvmi(v,m,i) for_wvmj(v,m,j)
#define for_wvmik(v,m,i,k) for_wvmi(v,m,i) for_wvmk(v,m,k)
#define for_wvmjk(v,m,j,k) for_wvmj(v,m,j) for_wvmk(v,m,k)

#endif

//----
// 3D  
//----
#define for_mijk(m,i,j,k) for_mi(m,i) for_mj(m,j) for_mk(m,k)

/* a */
#define for_amijk(m,i,j,k) for_ami(m,i) for_amj(m,j) for_amk(m,k)

/* v */
#define for_vmijk(v,m,i,j,k) for_vmi(v,m,i) for_vmj(v,m,j) for_vmk(v,m,k)

/* av */
#define for_avmijk(v,m,i,j,k) for_avmi(v,m,i) for_avmj(v,m,j) for_avmk(v,m,k)

#if 1
/* w */
#define for_wmijk(m,i,j,k) for_wmi(m,i) for_wmj(m,j) for_wmk(m,k)

/* wv */
#define for_wvmijk(v,m,i,j,k) for_wvmi(v,m,i) for_wvmj(v,m,j) for_wvmk(v,m,k)

#endif
#endif
