#include "vector.h"

/******************************************************************************/
Vector::Vector(const Domain & d) : aliav(false), vec() {
/*---------------------------------------------------------------------+
|  this constructor is used to decompose the vector over communicator  |
+---------------------------------------------------------------------*/
	
  dom = & d;
	
  /* assign this resolution to each vector component 
     and correct the values for staggered grid */
  for_m(m) {
    vec[m].n_x = ni();
    vec[m].n_y = nj();
    vec[m].n_z = nk();
  }
  vec[Comp::u()].n_x++;
  vec[Comp::v()].n_y++;
  vec[Comp::w()].n_z++;
  
  /* allocate all vector components */
  for_m(m) 
    vec[m].allocate( vec[m].n_x, vec[m].n_y, vec[m].n_z );
 
  /* correct offsets for each vector component */
  vec[Comp::u()].o_x++;
  vec[Comp::v()].o_y++;
  vec[Comp::w()].o_z++;

  /* correct offsets for each vector component */
#if 0
  vec[Comp::u()].s_X = 1;
  vec[Comp::u()].e_X = d.gi()-1;
  vec[Comp::u()].s_Y = 1;
  vec[Comp::u()].e_Y = d.gj()-2;
  vec[Comp::u()].s_Z = 1;
  vec[Comp::u()].e_Z = d.gk()-2;

  vec[Comp::v()].s_X = 1;
  vec[Comp::v()].e_X = d.gi()-2;
  vec[Comp::v()].s_Y = 1;
  vec[Comp::v()].e_Y = d.gj()-1;
  vec[Comp::v()].s_Z = 1;
  vec[Comp::v()].e_Z = d.gk()-2;

  vec[Comp::w()].s_X = 1;
  vec[Comp::w()].e_X = d.gi()-2;
  vec[Comp::w()].s_Y = 1;
  vec[Comp::w()].e_Y = d.gj()-2;
  vec[Comp::w()].s_Z = 1;
  vec[Comp::w()].e_Z = d.gk()-1;
#endif

  /* allocate memory for boundary conditions */
  for_m(m) 
    vec[m].bndcnd = new BndCnd( *dom );

  vec[Comp::u()].dom    = dom;
  vec[Comp::v()].dom    = dom;
  vec[Comp::w()].dom    = dom;

  coordinate();
}	

/******************************************************************************/
Vector::Vector(const Vector * v) : aliav(true), vec(v->vec) {
/*------------------------------------+
|  this constructor creates an alias  |
+------------------------------------*/

  for_m(m) 
    vec[m] = v->vec[m];

  for_m(m) {
    vec[m].n_x = v->vec[m].n_x; 
    vec[m].s_x = v->vec[m].s_x; 
    //vec[m].s_X = v->vec[m].s_X; 
    vec[m].e_x = v->vec[m].e_x;
    //vec[m].e_X = v->vec[m].e_X;
    vec[m].o_x = v->vec[m].o_x;
    vec[m].n_y = v->vec[m].n_y; 
    vec[m].s_y = v->vec[m].s_y; 
    //vec[m].s_Y = v->vec[m].s_Y; 
    vec[m].e_y = v->vec[m].e_y;
    //vec[m].e_Y = v->vec[m].e_Y;
    vec[m].o_y = v->vec[m].o_y;
    vec[m].n_z = v->vec[m].n_z; 
    vec[m].s_z = v->vec[m].s_z; 
    //vec[m].s_Z = v->vec[m].s_Z; 
    vec[m].e_z = v->vec[m].e_z;
    //vec[m].e_Z = v->vec[m].e_Z;
    vec[m].o_z = v->vec[m].o_z;

    pnt_xc [~m] = v->pnt_xc [~m]; 
    pnt_yc [~m] = v->pnt_yc [~m]; 
    pnt_zc [~m] = v->pnt_zc [~m]; 

    pnt_xn [~m] = v->pnt_xn [~m]; 
    pnt_yn [~m] = v->pnt_yn [~m]; 
    pnt_zn [~m] = v->pnt_zn [~m]; 

    pnt_xc_global [~m] = v->pnt_xc_global [~m];
    pnt_yc_global [~m] = v->pnt_yc_global [~m];
    pnt_zc_global [~m] = v->pnt_zc_global [~m];

    pnt_xn_global [~m] = v->pnt_xn_global [~m];
    pnt_yn_global [~m] = v->pnt_yn_global [~m];
    pnt_zn_global [~m] = v->pnt_zn_global [~m];

    pnt_dxc[~m] = v->pnt_dxc[~m]; 
    pnt_dyc[~m] = v->pnt_dyc[~m]; 
    pnt_dzc[~m] = v->pnt_dzc[~m]; 

    pnt_dxw[~m] = v->pnt_dxw[~m]; 
    pnt_dxe[~m] = v->pnt_dxe[~m]; 
    pnt_dys[~m] = v->pnt_dys[~m]; 
    pnt_dyn[~m] = v->pnt_dyn[~m]; 
    pnt_dzb[~m] = v->pnt_dzb[~m];
    pnt_dzt[~m] = v->pnt_dzt[~m];
  }

  dom    = v->dom;

  for_m(m)
    vec[m].bndcnd = v->vec[m].bndcnd;
}	

/******************************************************************************/
void Vector::allocate(const int ni, const int nj, const int nk) {
	
  /* assign this resolution to each vector component */
  for_m(m) {
    vec[m].n_x = ni;          
    vec[m].n_y = nj;       
    vec[m].n_z = nk;
  }
  vec[Comp::u()].n_x++;
  vec[Comp::v()].n_y++;
  vec[Comp::w()].n_z++;

  for_m(m) 
    vec[m].allocate( vec[m].n_x, vec[m].n_y, vec[m].n_z );

  /* correct offsets for each vector component */
  vec[Comp::u()].o_x++;
  vec[Comp::v()].o_y++;
  vec[Comp::w()].o_z++;
}	

/******************************************************************************/
void Vector::coordinate() {
/*-----------------------------------------------------+
|  sets the coordinate and cell dimension pointers up  |
+-----------------------------------------------------*/
	
  /* xc */
  pnt_xc[0]  = &Vector:: xc_staggered;
  pnt_xc[1]  = &Vector:: xc_nrm;
  pnt_xc[2]  = &Vector:: xc_nrm;
	
  /* yc */
  pnt_yc[0]  = &Vector:: yc_nrm;
  pnt_yc[1]  = &Vector:: yc_staggered;
  pnt_yc[2]  = &Vector:: yc_nrm;
	
  /* zc */
  pnt_zc[0]  = &Vector:: zc_nrm;
  pnt_zc[1]  = &Vector:: zc_nrm;
  pnt_zc[2]  = &Vector:: zc_staggered;
	
  /* xn */
  pnt_xn[0]  = &Vector:: xn_staggered;
  pnt_xn[1]  = &Vector:: xn_nrm;
  pnt_xn[2]  = &Vector:: xn_nrm;
	
  /* yn */
  pnt_yn[0]  = &Vector:: yn_nrm;
  pnt_yn[1]  = &Vector:: yn_staggered;
  pnt_yn[2]  = &Vector:: yn_nrm;
	
  /* zn */
  pnt_zn[0]  = &Vector:: zn_nrm;
  pnt_zn[1]  = &Vector:: zn_nrm;
  pnt_zn[2]  = &Vector:: zn_staggered;

  /* xc_global */
  pnt_xc_global[0]  = &Vector:: xc_staggered_global;
  pnt_xc_global[1]  = &Vector:: xc_nrm_global;
  pnt_xc_global[2]  = &Vector:: xc_nrm_global;

  /* yc_global */
  pnt_yc_global[0]  = &Vector:: yc_nrm_global;
  pnt_yc_global[1]  = &Vector:: yc_staggered_global;
  pnt_yc_global[2]  = &Vector:: yc_nrm_global;

  /* zc_global */
  pnt_zc_global[0]  = &Vector:: zc_nrm_global;
  pnt_zc_global[1]  = &Vector:: zc_nrm_global;
  pnt_zc_global[2]  = &Vector:: zc_staggered_global;

  /* xn_global */
  pnt_xn_global[0]  = &Vector:: xn_staggered_global;
  pnt_xn_global[1]  = &Vector:: xn_nrm_global;
  pnt_xn_global[2]  = &Vector:: xn_nrm_global;

  /* yn_global */
  pnt_yn_global[0]  = &Vector:: yn_nrm_global;
  pnt_yn_global[1]  = &Vector:: yn_staggered_global;
  pnt_yn_global[2]  = &Vector:: yn_nrm_global;

  /* zn_global */
  pnt_zn_global[0]  = &Vector:: zn_nrm_global;
  pnt_zn_global[1]  = &Vector:: zn_nrm_global;
  pnt_zn_global[2]  = &Vector:: zn_staggered_global;
	
  /* dxc */
  pnt_dxc[0] = &Vector::dxc_staggered;
  pnt_dxc[1] = &Vector::dxc_nrm;
  pnt_dxc[2] = &Vector::dxc_nrm;
	
  pnt_dxe[0] = &Vector::dxe_staggered;
  pnt_dxe[1] = &Vector::dxe_nrm;
  pnt_dxe[2] = &Vector::dxe_nrm;
	
  pnt_dxw[0] = &Vector::dxw_staggered;
  pnt_dxw[1] = &Vector::dxw_nrm;
  pnt_dxw[2] = &Vector::dxw_nrm;
	
  /* dyc */
  pnt_dyc[0] = &Vector::dyc_nrm;
  pnt_dyc[1] = &Vector::dyc_staggered;
  pnt_dyc[2] = &Vector::dyc_nrm;
	
  pnt_dys[0] = &Vector::dys_nrm;
  pnt_dys[1] = &Vector::dys_staggered;
  pnt_dys[2] = &Vector::dys_nrm;
	
  pnt_dyn[0] = &Vector::dyn_nrm;
  pnt_dyn[1] = &Vector::dyn_staggered;
  pnt_dyn[2] = &Vector::dyn_nrm;
	
  /* dzc */
  pnt_dzc[0] = &Vector::dzc_nrm;
  pnt_dzc[1] = &Vector::dzc_nrm;
  pnt_dzc[2] = &Vector::dzc_staggered;
	
  pnt_dzb[0] = &Vector::dzb_nrm;
  pnt_dzb[1] = &Vector::dzb_nrm;
  pnt_dzb[2] = &Vector::dzb_staggered;
	
  pnt_dzt[0] = &Vector::dzt_nrm;
  pnt_dzt[1] = &Vector::dzt_nrm;
  pnt_dzt[2] = &Vector::dzt_staggered;
	
}	

/******************************************************************************/
Vector::~Vector() {

  if(!aliav)
    for_m(m)
      vec[m].deallocate();
}	

/*-----------------------------------------------------------------------------+
 '$Id: vector.cpp,v 1.20 2016/03/15 15:36:34 sato Exp $'/
+-----------------------------------------------------------------------------*/
