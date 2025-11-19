#ifndef BROWSING_LAGRANGIAN_H
#define BROWSING_LAGRANGIAN_H

//-------------------
// through particles 
//-------------------
#define for_p(p)    for(int p=0; p<size(); p++)   /* use inside the class */

#define for_dp(d,p) for(int p=0; p<d.size(); p++) /* use outside the class */

//----
// 1D
//----
#define for_pi(p,i) for(int i=particles[p].si(); i<=particles[p].ei(); i++)
#define for_pj(p,j) for(int j=particles[p].sj(); j<=particles[p].ej(); j++)
#define for_pk(p,k) for(int k=particles[p].sk(); k<=particles[p].ek(); k++)

#define for_dpi(d,p,i) for(int i=d[p].si(); i<=d[p].ei(); i++)
#define for_dpj(d,p,j) for(int j=d[p].sj(); j<=d[p].ej(); j++)
#define for_dpk(d,p,k) for(int k=d[p].sk(); k<=d[p].ek(); k++)

//----
// 3D
//----
#define for_pijk(p,i,j,k) for_pi(p,i) for_pj(p,j) for_pk(p,k)

#define for_dpijk(d,p,i,j,k) for_dpi(d,p,i) for_dpj(d,p,j) for_dpk(d,p,k)

#endif

/*-----------------------------------------------------------------------------+
 '$Id: lagrangian_browsing.h,v 1.1 2018/06/15 11:39:11 MinZhu Exp $'/ 
+-----------------------------------------------------------------------------*/
