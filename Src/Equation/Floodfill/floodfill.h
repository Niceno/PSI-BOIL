#ifndef FLOODFILL_H
#define FLOODFILL_H

#include <vector>
#include "index_ijk.h"
#include "../../Field/Vector/vector.h"
#include "../../Plot/plot.h"
#include "../../SimulationTime/simulation_time.h"
#include "region.h"


/////////////////
//             //
//  Floodfill  //
//             //
/////////////////
class Floodfill { 

  public:
    Floodfill(Scalar & colorf, 
              Scalar & srid, 
              const Vector * const p_uvw,
              const Times & t);
    ~Floodfill();

    /* main function identify all closed regions in domain */
    void identify_regions();

    void save(const char * nm, const int it);
    void load(const char * nm, const int it);

    void set_out_freq(const int i){
      out_rgn_info_freq=i;
      boil::oout<<"Floodfill:out_freq= "<<i<<"\n";
    }
    int  get_out_freq(){return out_rgn_info_freq;}

  private:
    std::vector<Index_ijk> new_seed_stack, same_fill_stack, track_index_stack;
    std::vector<int> rgnid_match_arrayi, rgnid_match_arrayj,
                     rgnid_match_array_gatheredi, rgnid_match_array_gatheredj;
    std::vector<int> lookup_rgnid_array, smlookup_oldid, smlookupvchk, oldid_num_newrgns,
                     lookup_volcells, smlookup_cvol, flookup_cvol;
    std::vector< std::vector<int> > track_rid_matrix;
    std::vector<real> flookup_x, flookup_y, flookup_z,
                      flookup_u, flookup_v, flookup_w;

    std::vector<int> v_idavail;
    int * pt_idavail; 
 
    const Times * time;
    real oldtime; //for outputting dx(com)/dt over longer dt

    int m_rgn_n;
    int m_rgn_p;

    int old_neg_smrgnid;
    int old_pos_smrgnid;

    /* color function criterion to check for interface */
    const real c_criterion = 0.5;
    
    const real m_dV;
    const real m_dxc;

    int prgn_inc;
    int nrgn_inc;

    /* for labeling debug plotting */
    int pltinc;

    Scalar & c;
    const Vector * const uvw;
    Scalar & rgnid;
    Scalar rgnid_old;

    /* for parallel */
    int * num_rgn_pairs_gth;
    int * displace_index;

    void scanline();

    void compare_rgn(const int & rgnin, const int & rgnout);

    void face_scan_rgn();

    void get_region_info(int ptmpcalcID, int ntmpcalcID);
    void out_region_info();

    /* have an array of regions */
    std::vector<Region> m_vectrgns;
 
    int * pt_flookup_cvol;
    real * pt_flookup_x;
    real * pt_flookup_y;
    real * pt_flookup_z;
    real * pt_flookup_u;
    real * pt_flookup_v;
    real * pt_flookup_w;

    Region & getregion(int rid);  //was posregion or negregion
    
    /* rgn erased after hidden for more than 30 time steps (why 30, good guess) */
    int mng_hidden_rgns();  
    
    int hidden_rgnid(real ix, real iy, real iz); //was hid_prgn_xyz, hid_nrgn_xyz

    void  erase_region(int rid);

    int get_new_pos_id();
    int get_new_neg_id();

    int out_rgn_info_freq;
    
    std::ofstream outrgn;

};

#endif
