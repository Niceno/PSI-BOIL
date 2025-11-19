#include "vectorbool.h"

/******************************************************************************/
VectorBool::VectorBool(const Domain & d) : aliav(false), vec() {
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

  /* allocate memory for boundary conditions */
  for_m(m) 
    vec[m].bndcnd = new BndCnd( *dom );

  vec[Comp::u()].dom    = dom;
  vec[Comp::v()].dom    = dom;
  vec[Comp::w()].dom    = dom;

  coordinate();
}	
#if 1
/******************************************************************************/
VectorBool::VectorBool(const VectorBool * v) : aliav(true), vec(v->vec) {
/*------------------------------------+
|  this constructor creates an alias  |
+------------------------------------*/

  for_m(m) 
    vec[m] = v->vec[m];

  for_m(m) {
    vec[m].n_x = v->vec[m].n_x; 
    vec[m].s_x = v->vec[m].s_x; 
    vec[m].e_x = v->vec[m].e_x;
    vec[m].o_x = v->vec[m].o_x;
    vec[m].n_y = v->vec[m].n_y; 
    vec[m].s_y = v->vec[m].s_y; 
    vec[m].e_y = v->vec[m].e_y;
    vec[m].o_y = v->vec[m].o_y;
    vec[m].n_z = v->vec[m].n_z; 
    vec[m].s_z = v->vec[m].s_z; 
    vec[m].e_z = v->vec[m].e_z;
    vec[m].o_z = v->vec[m].o_z;

    pnt_xc [~m] = v->pnt_xc [~m]; 
    pnt_yc [~m] = v->pnt_yc [~m]; 
    pnt_zc [~m] = v->pnt_zc [~m]; 

    pnt_xn [~m] = v->pnt_xn [~m]; 
    pnt_yn [~m] = v->pnt_yn [~m]; 
    pnt_zn [~m] = v->pnt_zn [~m]; 

    pnt_dxc[~m] = v->pnt_dxc[~m]; 
    pnt_dyc[~m] = v->pnt_dyc[~m]; 
    pnt_dzc[~m] = v->pnt_dzc[~m]; 

    pnt_dxw[~m] = v->pnt_dxw[~m]; 
    pnt_dxe[~m] = v->pnt_dxe[~m]; 
    pnt_dys[~m] = v->pnt_dys[~m]; 
    pnt_dyn[~m] = v->pnt_dyn[~m]; 
    pnt_dzb[~m] = v->pnt_dzb[~m];
    pnt_dzt[~m] = v->pnt_dzt[~m];

    pnt_dSx[~m] = v->pnt_dSx[~m];
    pnt_dSy[~m] = v->pnt_dSy[~m];
    pnt_dSz[~m] = v->pnt_dSz[~m];

    pnt_dSx_dir[~m] = v->pnt_dSx_dir[~m];
    pnt_dSy_dir[~m] = v->pnt_dSy_dir[~m];
    pnt_dSz_dir[~m] = v->pnt_dSz_dir[~m];

    pnt_dV[~m] = v->pnt_dV[~m];
  }

  dom    = v->dom;

  for_m(m)
    vec[m].bndcnd = v->vec[m].bndcnd;
}	
#endif
/******************************************************************************/
void VectorBool::allocate(const int ni, const int nj, const int nk) {
	
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
void VectorBool::coordinate() {
/*-----------------------------------------------------+
|  sets the coordinate and cell dimension pointers up  |
+-----------------------------------------------------*/
	
  /* xc */
  pnt_xc[0]  = &VectorBool:: xc_staggered;
  pnt_xc[1]  = &VectorBool:: xc_nrm;
  pnt_xc[2]  = &VectorBool:: xc_nrm;
	
  /* yc */
  pnt_yc[0]  = &VectorBool:: yc_nrm;
  pnt_yc[1]  = &VectorBool:: yc_staggered;
  pnt_yc[2]  = &VectorBool:: yc_nrm;
	
  /* zc */
  pnt_zc[0]  = &VectorBool:: zc_nrm;
  pnt_zc[1]  = &VectorBool:: zc_nrm;
  pnt_zc[2]  = &VectorBool:: zc_staggered;
	
  /* xn */
  pnt_xn[0]  = &VectorBool:: xn_staggered;
  pnt_xn[1]  = &VectorBool:: xn_nrm;
  pnt_xn[2]  = &VectorBool:: xn_nrm;
	
  /* yn */
  pnt_yn[0]  = &VectorBool:: yn_nrm;
  pnt_yn[1]  = &VectorBool:: yn_staggered;
  pnt_yn[2]  = &VectorBool:: yn_nrm;
	
  /* zn */
  pnt_zn[0]  = &VectorBool:: zn_nrm;
  pnt_zn[1]  = &VectorBool:: zn_nrm;
  pnt_zn[2]  = &VectorBool:: zn_staggered;
	
  /* dxc */
  pnt_dxc[0] = &VectorBool::dxc_staggered;
  pnt_dxc[1] = &VectorBool::dxc_nrm;
  pnt_dxc[2] = &VectorBool::dxc_nrm;
	
  pnt_dxe[0] = &VectorBool::dxe_staggered;
  pnt_dxe[1] = &VectorBool::dxe_nrm;
  pnt_dxe[2] = &VectorBool::dxe_nrm;
	
  pnt_dxw[0] = &VectorBool::dxw_staggered;
  pnt_dxw[1] = &VectorBool::dxw_nrm;
  pnt_dxw[2] = &VectorBool::dxw_nrm;
	
  /* dyc */
  pnt_dyc[0] = &VectorBool::dyc_nrm;
  pnt_dyc[1] = &VectorBool::dyc_staggered;
  pnt_dyc[2] = &VectorBool::dyc_nrm;
	
  pnt_dys[0] = &VectorBool::dys_nrm;
  pnt_dys[1] = &VectorBool::dys_staggered;
  pnt_dys[2] = &VectorBool::dys_nrm;
	
  pnt_dyn[0] = &VectorBool::dyn_nrm;
  pnt_dyn[1] = &VectorBool::dyn_staggered;
  pnt_dyn[2] = &VectorBool::dyn_nrm;
	
  /* dzc */
  pnt_dzc[0] = &VectorBool::dzc_nrm;
  pnt_dzc[1] = &VectorBool::dzc_nrm;
  pnt_dzc[2] = &VectorBool::dzc_staggered;
	
  pnt_dzb[0] = &VectorBool::dzb_nrm;
  pnt_dzb[1] = &VectorBool::dzb_nrm;
  pnt_dzb[2] = &VectorBool::dzb_staggered;
	
  pnt_dzt[0] = &VectorBool::dzt_nrm;
  pnt_dzt[1] = &VectorBool::dzt_nrm;
  pnt_dzt[2] = &VectorBool::dzt_staggered;
	
  /* dS */
  pnt_dSx[0] = &VectorBool::dSx_xstaggered;
  pnt_dSx[1] = &VectorBool::dSx_ystaggered;
  pnt_dSx[2] = &VectorBool::dSx_zstaggered;

  pnt_dSy[0] = &VectorBool::dSy_xstaggered;
  pnt_dSy[1] = &VectorBool::dSy_ystaggered;
  pnt_dSy[2] = &VectorBool::dSy_zstaggered;

  pnt_dSz[0] = &VectorBool::dSz_xstaggered;
  pnt_dSz[1] = &VectorBool::dSz_ystaggered;
  pnt_dSz[2] = &VectorBool::dSz_zstaggered;

  /* dS dir */
  pnt_dSx_dir[0] = &VectorBool::dSx_dir_xstaggered;
  pnt_dSx_dir[1] = &VectorBool::dSx_dir_ystaggered;
  pnt_dSx_dir[2] = &VectorBool::dSx_dir_zstaggered;

  pnt_dSy_dir[0] = &VectorBool::dSy_dir_xstaggered;
  pnt_dSy_dir[1] = &VectorBool::dSy_dir_ystaggered;
  pnt_dSy_dir[2] = &VectorBool::dSy_dir_zstaggered;

  pnt_dSz_dir[0] = &VectorBool::dSz_dir_xstaggered;
  pnt_dSz_dir[1] = &VectorBool::dSz_dir_ystaggered;
  pnt_dSz_dir[2] = &VectorBool::dSz_dir_zstaggered;

  /* dV */
  pnt_dV[0] = &VectorBool::dV_xstaggered;
  pnt_dV[1] = &VectorBool::dV_ystaggered;
  pnt_dV[2] = &VectorBool::dV_zstaggered;

  return;
}	

/******************************************************************************/
VectorBool::~VectorBool() {

  if(!aliav)
    for_m(m)
      vec[m].deallocate();
}
