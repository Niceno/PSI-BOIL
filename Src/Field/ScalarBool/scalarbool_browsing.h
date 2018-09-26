#ifndef BROWSING_SCALAR_H
#define BROWSING_SCALAR_H

//----
// 1D
//----
#define for_i(i) for(int i=si(); i<=ei(); i++)
#define for_j(j) for(int j=sj(); j<=ej(); j++)
#define for_k(k) for(int k=sk(); k<=ek(); k++)

/* a */
#define for_ai(i) for(int i=0; i<ni(); i++)
#define for_aj(j) for(int j=0; j<nj(); j++)
#define for_ak(k) for(int k=0; k<nk(); k++)

/* av */
#define for_avi(v,i) for(int i=0; i<v.ni(); i++)
#define for_avj(v,j) for(int j=0; j<v.nj(); j++)
#define for_avk(v,k) for(int k=0; k<v.nk(); k++)

/* v */
#define for_vi(v,i) for(int i=v.si(); i<=v.ei(); i++)
#define for_vj(v,j) for(int j=v.sj(); j<=v.ej(); j++)
#define for_vk(v,k) for(int k=v.sk(); k<=v.ek(); k++)

//----
// 2D
//----
#define for_ij(i,j) for_i(i) for_j(j)
#define for_ik(i,k) for_i(i) for_k(k)
#define for_jk(j,k) for_j(j) for_k(k)

/* a */
#define for_aij(i,j) for_ai(i) for_aj(j)
#define for_aik(i,k) for_ai(i) for_ak(k)
#define for_ajk(j,k) for_aj(j) for_ak(k)

/* v */
#define for_avij(v,i,j) for_avi(v,i) for_avj(v,j)
#define for_avik(v,i,k) for_avi(v,i) for_avk(v,k)
#define for_avjk(v,j,k) for_avj(v,j) for_avk(v,k)

/* v */
#define for_vij(v,i,j) for_vi(v,i) for_vj(v,j)
#define for_vik(v,i,k) for_vi(v,i) for_vk(v,k)
#define for_vjk(v,j,k) for_vj(v,j) for_vk(v,k)

//----
// 3D
//----
#define for_ijk(i,j,k) for_i(i) for_j(j) for_k(k)

/* a */
#define for_aijk(i,j,k) for_ai(i) for_aj(j) for_ak(k)

/* av */
#define for_avijk(v,i,j,k) for_avi(v,i) for_avj(v,j) for_avk(v,k)

/* v */
#define for_vijk(v,i,j,k) for_vi(v,i) for_vj(v,j) for_vk(v,k)

#endif

/*-----------------------------------------------------------------------------+
 '$Id: scalarbool_browsing.h,v 1.1 2014/02/04 08:16:57 sato Exp $'/
+-----------------------------------------------------------------------------*/
