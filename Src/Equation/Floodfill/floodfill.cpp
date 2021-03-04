#include "floodfill.h"

#define BIG_BUFFER_SIZE   (262144)
#define SMALL_BUFFER_SIZE (BIG_BUFFER_SIZE/8)

/******************************************************************************
  This is a non-trivial constructor. 
*******************************************************************************/
Floodfill::Floodfill(Scalar & colorf, 
                     Scalar & srid, 
                     const Vector * const p_uvw,
                     const Times  & t ) :
                     /* constructor initialization list */
                     c(colorf), 
                     rgnid(srid), 
                     uvw(p_uvw),
                     time ( & t),
                     //rgnid_old(*(srid.domain())) {
                     rgnid_old(*(srid.domain())),
                     //m_dV(srid.dV(srid.si()+1,srid.sj()+1,srid.sk()+1)),
                     //m_dxc(srid.dxc(srid.si()+1)) {
                     m_dV(srid.domain()->dV_min()),
                     m_dxc(srid.domain()->dxyz_min()) {

  new_seed_stack.reserve (BIG_BUFFER_SIZE);  /* sizes may be too small/big? */
  same_fill_stack.reserve(BIG_BUFFER_SIZE);
  track_index_stack.reserve(BIG_BUFFER_SIZE);
  rgnid_match_arrayi.reserve(BIG_BUFFER_SIZE);
  rgnid_match_arrayj.reserve(BIG_BUFFER_SIZE);
  rgnid_match_array_gatheredi.reserve(boil::cart.nproc()*SMALL_BUFFER_SIZE);
  rgnid_match_array_gatheredj.reserve(boil::cart.nproc()*SMALL_BUFFER_SIZE);
  lookup_rgnid_array.reserve(BIG_BUFFER_SIZE);
  lookup_volcells.reserve(BIG_BUFFER_SIZE);
  flookup_x.reserve(BIG_BUFFER_SIZE);
  flookup_y.reserve(BIG_BUFFER_SIZE);
  flookup_z.reserve(BIG_BUFFER_SIZE);
  flookup_u.reserve(BIG_BUFFER_SIZE);
  flookup_v.reserve(BIG_BUFFER_SIZE);
  flookup_w.reserve(BIG_BUFFER_SIZE);
  flookup_cvol.reserve(BIG_BUFFER_SIZE);
  smlookup_oldid.reserve(BIG_BUFFER_SIZE);
  smlookupvchk.reserve(BIG_BUFFER_SIZE);
  oldid_num_newrgns.reserve(BIG_BUFFER_SIZE);
  smlookup_cvol.reserve(BIG_BUFFER_SIZE);

  /* idavail vector initially contains 1 element id=0 and false for not available
     because region id 0 should never be used */
  v_idavail.reserve(SMALL_BUFFER_SIZE);
  v_idavail.resize(1, false); 
  pt_idavail = v_idavail.data();

  pt_flookup_cvol = 0;
  pt_flookup_x   = 0;
  pt_flookup_y   = 0;
  pt_flookup_z   = 0;
  pt_flookup_u   = 0;
  pt_flookup_v   = 0;
  pt_flookup_w   = 0;

  rgnid_old = srid.shape();
  
  prgn_inc = 0; //laff important not to redeclare these!!!, just assign!
  nrgn_inc = 0;
  pltinc = 0; //for debug plotting

  out_rgn_info_freq=50;

  if (!boil::cart.iam()){
    outrgn.open("tracked_regions.txt",std::ios::app);
    boil::oout<<"floodfill: open file 'tracked_regions.txt'\n";
  }

  old_neg_smrgnid = 0;
  old_pos_smrgnid = 0;
  oldtime = 0.0;

  num_rgn_pairs_gth = new int[boil::cart.nproc()](); /* zero initialize IAW   */
  displace_index = new int[boil::cart.nproc()]();    /* ISO C++ Section 8.5/5 */

  boil::oout<<"Floodfill: output will be stored in 'tracked_regions.txt'.\n";
}
 

/******************************************************************************/
Floodfill::~Floodfill(){
  /*-------------------------------------------------------------+
  |  destructor to deallocate memory and reset pointers to null  |
  |  for dynamic memory management                               |
  +-------------------------------------------------------------*/
  delete[] num_rgn_pairs_gth; /* deallocate memory */
  num_rgn_pairs_gth = 0;      /* reset pointer to null */
  delete[] displace_index;
  displace_index = 0;    
  if (!boil::cart.iam())
    outrgn.close();
}

/******************************************************************************/
Region & Floodfill::getregion(int rid) {
  for (int i=0; i<m_vectrgns.size(); i++) {
    if (m_vectrgns[i].id() == rid) return m_vectrgns[i];
  }
  m_vectrgns.push_back(Region(rid));
  return m_vectrgns.back();
}

/******************************************************************************/
int Floodfill::mng_hidden_rgns() {
  /*------------------------------------------------------------------+
  |  rgn erased after hidden for more than 100(why?) time steps       |
  +------------------------------------------------------------------*/
  const int hidden_limit = 100;
  for (int i=0; i<m_vectrgns.size(); i++) {
    if (m_vectrgns[i].hiding() == true) {
      if (m_vectrgns[i].inc_hiding() > hidden_limit) {
        Region & r = m_vectrgns[i];
        if(!boil::cart.iam()) {
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<r.id()<<" cvol "
          <<r.cellvol()<<" xyz "<<r.x()<<" "<<r.y()<<" "<<r.z()
          <<" in hiding is being erased"<<std::endl;
        }
        pt_idavail[m_vectrgns[i].id()] = true;
        m_vectrgns.erase(m_vectrgns.begin()+i);
      }
    }
  }
  return 0;
}

/******************************************************************************/
int Floodfill::hidden_rgnid(real ix, real iy, real iz) { 
  for (int i=0; i<m_vectrgns.size(); i++) {
    if (m_vectrgns[i].hiding() == true) {
      real idx = (m_vectrgns[i].x()-ix);
      real idy = (m_vectrgns[i].y()-iy);
      real idz = (m_vectrgns[i].z()-iz);
      real dist = sqrt(idx*idx + idy*idy + idz*idz);
      if (dist < 1.5*m_dxc) {
        if (!boil::cart.iam()) {
          outrgn<<"t "<<time->current_time()<<" tsteps in hiding "
                <<m_vectrgns[i].get_tsteps_hidden()<<" dist traveled in hiding "
                <<dist<<" < 1.5dx "<<1.5*m_dxc<<std::endl;
        }
        m_vectrgns[i].hiding(false);
        return m_vectrgns[i].id();
      }
    }
  }
  return 0; 
}

/******************************************************************************/
void Floodfill::erase_region(int rid){
  /*------------------------------------------------------------------+
  |  erase region from vector of regions with id=rid                  |
  +------------------------------------------------------------------*/
  for (int i=0; i<m_vectrgns.size(); i++) {
    if (m_vectrgns[i].id() == rid) {
      pt_idavail[rid] = true;
      m_vectrgns.erase(m_vectrgns.begin()+i);
    }
  }
}

/******************************************************************************/
int Floodfill::get_new_pos_id() {
  for (int i=1; i<=prgn_inc; i++) { 
    if (pt_idavail[i] == true) {
      pt_idavail[i] = false;
      return i; 
    }
  }
  v_idavail.push_back(false);
  prgn_inc++;
  return prgn_inc;
}

/******************************************************************************/
int Floodfill::get_new_neg_id() {
  for (int i=-1; i>=nrgn_inc; i--) { 
    if (pt_idavail[i] == true) {
      pt_idavail[i] = false;
      return i;
    }
  }
  v_idavail.insert(v_idavail.begin(), false);
  nrgn_inc--;
  pt_idavail = v_idavail.data() - nrgn_inc;
  return nrgn_inc;
}
