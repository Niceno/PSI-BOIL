#include <iomanip>
#include "floodfill.h"
#include "../../Plot/plot.h"
//#define VERBOSE

/***************************************************************************//**
*  This function should be called from the main.cpp and calls other functions
*  to handle the floodfill and parallel tasks. 
*******************************************************************************/
void Floodfill::identify_regions() {

  boil::timer.start("id_rgns_mn");

  /* for debug plotting when merging or splitting */
  bool need2plot = false;

  /* if region is small may not be tracked properly ( could not
     be detected and then be redetected if smeared color function drops below
     0.5 cutoff */
  const int size_smallrgn = 3;

  pt_idavail = v_idavail.data() - nrgn_inc;


  for_vijk(rgnid,i,j,k) {
    rgnid_old[i][j][k] = rgnid[i][j][k];
    rgnid[i][j][k] = 0;
  }
  rgnid_old.exchange();
  rgnid.exchange();

  old_pos_smrgnid = prgn_inc;
  old_neg_smrgnid = nrgn_inc;

  
  if (!new_seed_stack.empty()) {
    boil::aout<<"EXIT: new_seed_stack not empty in FF_identify_regions()"
      <<__LINE__<< boil::endl;
    exit(1);
  }

  new_seed_stack.push_back(Index_ijk(rgnid.si()+1,rgnid.sj()+1,rgnid.sk()+1));
  
  m_rgn_n = 0;
  m_rgn_p = 0;
  int max_rgnid = 0;
  int min_rgnid = 0;

  /*--------------------------------------+
  |  Scan/Floodfill each region in domain |
  |  and give each region an ID           |
  |  that is (-) if c <= c_criterion and  |
  |  (+) if c > c_criterion               |
  +--------------------------------------*/
  while (!new_seed_stack.empty()) { 
    boil::timer.start("scanline");
    scanline();  
    boil::timer.stop("scanline");
  }

  int rgnid_ct_max = std::max(std::abs(m_rgn_n),std::abs(m_rgn_p));
  
  boil::cart.max_int_n(&rgnid_ct_max, 1);
  
  for_vijk(rgnid,i,j,k) {
    rgnid[i][j][k] = ((rgnid[i][j][k]>0) - (rgnid[i][j][k]<0))
                       * rgnid_ct_max * boil::cart.iam() 
                       + rgnid[i][j][k];
  }
  
  rgnid.exchange_all();
  
  /* empty these arrays */
  while(!rgnid_match_arrayi.empty()){ 
    rgnid_match_arrayi.pop_back();
  }
  while(!rgnid_match_arrayj.empty()){ 
    rgnid_match_arrayj.pop_back();
  }
  while(!rgnid_match_array_gatheredi.empty()) { 
    rgnid_match_array_gatheredi.pop_back(); 
  }
  while(!rgnid_match_array_gatheredj.empty()) { 
    rgnid_match_array_gatheredj.pop_back(); 
  }
  boil::timer.start("id_facescan");
  face_scan_rgn();
  boil::timer.stop("id_facescan");
  
  std::fill(num_rgn_pairs_gth, num_rgn_pairs_gth + boil::cart.nproc(), 0);

  int rgnid_match_array_size = rgnid_match_arrayi.size();
  num_rgn_pairs_gth[boil::cart.iam()] = rgnid_match_array_size; 
  boil::cart.sum_int_n(num_rgn_pairs_gth, boil::cart.nproc());
  
  int sum_displc = 0;
  for (int i=0; i<boil::cart.nproc(); i++){
    displace_index[i] = sum_displc;
    sum_displc += num_rgn_pairs_gth[i];
  }

  rgnid_match_array_gatheredi.resize(sum_displc);
  rgnid_match_array_gatheredj.resize(sum_displc);

  for (int i=0; i<rgnid_match_array_size; i++) {
    int j = i+displace_index[boil::cart.iam()];
    rgnid_match_array_gatheredi[j]=rgnid_match_arrayi[i];
    rgnid_match_array_gatheredj[j]=rgnid_match_arrayj[i];
  }
  boil::cart.sum_int_n(&rgnid_match_array_gatheredi[0], sum_displc);
  boil::cart.sum_int_n(&rgnid_match_array_gatheredj[0], sum_displc);

  max_rgnid = rgnid_ct_max * (boil::cart.nproc()-1) 
                     + rgnid_ct_max;
  min_rgnid = -max_rgnid;

  /* for 1 proc execution */
  if (boil::cart.nproc() <= 1) {
    max_rgnid = m_rgn_p;
    min_rgnid = m_rgn_n;
  }

  int size_lookups = max_rgnid - min_rgnid + 1;
  lookup_rgnid_array.resize(size_lookups);
  lookup_volcells.resize(size_lookups);

  for (int i=0; i<size_lookups; i++) {
    lookup_rgnid_array[i] = min_rgnid + i; /* initialize with rgnid */
    lookup_volcells[i] = 0;
  }

  int * pt_lookup_rgnid_array = lookup_rgnid_array.data();
  int * pt_lookup_volcells = lookup_volcells.data();
  pt_lookup_rgnid_array -= min_rgnid;
  pt_lookup_volcells    -= min_rgnid;

  /* quicksort algorithm performed so when a region id is looked up in the
     lookup_rgnid_array_gathered a reduced shared region ID is returned */
  boil::timer.start("id_quicksort");
  for (int i=0; i<rgnid_match_array_gatheredi.size(); i++) {
    int p = rgnid_match_array_gatheredi[i];
    int q = rgnid_match_array_gatheredj[i];
    int t = pt_lookup_rgnid_array[p];
    if (t == q) continue;
    for (int j=min_rgnid; j<=max_rgnid; j++) {
      if (pt_lookup_rgnid_array[j] == t)
        pt_lookup_rgnid_array[j] = pt_lookup_rgnid_array[q];
    }
  }
  boil::timer.stop("id_quicksort");

  for_vijk(rgnid,i,j,k) {
    int old_rgnid = int(rgnid[i][j][k]);
    rgnid[i][j][k] = pt_lookup_rgnid_array[old_rgnid];
  }
  
  for_vijk(rgnid,i,j,k) {  
    int rid = int(rgnid[i][j][k]);
    pt_lookup_volcells[rid]++;
  }
  boil::cart.sum_int_n(&lookup_volcells[0], size_lookups);

  int neg_smrgnid = 0;
  int pos_smrgnid = 0;
  for (int i=min_rgnid; i<=max_rgnid; i++) {
    if (pt_lookup_volcells[i] == 0) {
      pt_lookup_rgnid_array[i] = 0;
    } 
    else if (i < 0) {
      neg_smrgnid--;
      pt_lookup_rgnid_array[i] = neg_smrgnid;
    }
    else {
      pos_smrgnid++;
      pt_lookup_rgnid_array[i] = pos_smrgnid;
    }
  }
  
  for_vijk(rgnid,i,j,k) {
    int old_rgnid = int(rgnid[i][j][k]);
    rgnid[i][j][k] = pt_lookup_rgnid_array[old_rgnid];
  }


  int new_size_lookups = pos_smrgnid - neg_smrgnid + 1;

  smlookup_oldid.resize(new_size_lookups);
  smlookupvchk.resize(new_size_lookups);
  smlookup_cvol.resize(new_size_lookups);
  int * pt_smlookup_oldid = smlookup_oldid.data(); //input small tmpcalcId and get oldID or 0
  int * pt_smlookupvchk = smlookupvchk.data();
  int * pt_smlookup_cvol = smlookup_cvol.data();
  pt_smlookup_oldid -= neg_smrgnid;
  pt_smlookupvchk -= neg_smrgnid;
  pt_smlookup_cvol -= neg_smrgnid;

  for (int i=neg_smrgnid; i<=pos_smrgnid; i++) {
    pt_smlookup_oldid[i]=0;
    pt_smlookupvchk[i]=0;
    pt_smlookup_cvol[i]=0;
  }

  for (int i=min_rgnid; i<=max_rgnid; i++) {
    if (pt_lookup_rgnid_array[i]==0) {}//continue; //do nothing
    else {
      int smi = pt_lookup_rgnid_array[i];
      int voli = pt_lookup_volcells[i];
      pt_smlookup_cvol[smi] = voli;
    }
  }

  boil::timer.start("id_make_tr_matrix");
  int old_size_lookups = old_pos_smrgnid - old_neg_smrgnid + 1;
  track_rid_matrix.resize(new_size_lookups);
  for (int i=0; i<new_size_lookups; i++) {
    track_rid_matrix[i].resize(old_size_lookups);
  }
  
  /* initialize to 0 */
  for (int i=0; i<new_size_lookups; i++) {
    for (int j=0; j<old_size_lookups; j++) {
      track_rid_matrix[i][j] = 0;
    }
  }

  for_vijk(rgnid,i,j,k) {
    int rid = int(rgnid[i][j][k]) - neg_smrgnid;
    int rid_old = int(rgnid_old[i][j][k]) - old_neg_smrgnid;
    track_rid_matrix[rid][rid_old]++;
  }

  for (int i=0; i<new_size_lookups; i++) {
    boil::cart.sum_int_n(&track_rid_matrix[i][0], old_size_lookups);
  }
  boil::timer.stop("id_make_tr_matrix");

  get_region_info(pos_smrgnid, neg_smrgnid);

  oldid_num_newrgns.resize(old_size_lookups);
  int * pt_oldid_num_newrgns = oldid_num_newrgns.data();
  pt_oldid_num_newrgns -= old_neg_smrgnid;
  for (int i=old_neg_smrgnid; i<=old_pos_smrgnid; i++) {
    pt_oldid_num_newrgns[i]=0;
  }

  /*------------------------------------------------------------------------------+
  |  update new_rgnids (tmpcalcIDs) with largest matching old_rgnids              | 
  |  (0 otherwise, if smaller id from a split, or no matching old id)             |
  |  and  store number of new regions (tmpcaclIDs) originating (splitting?) from  |
  |  a given old region in pt_oldid_num_newrgns[]                                 |
  |  pt_oldid_num_newrgns[oldID]->number of new IDs from oldID                    |
  +------------------------------------------------------------------------------*/
  boil::timer.start("id_new->old");
  for (int j=(old_size_lookups-old_pos_smrgnid); j<old_size_lookups; j++) {
    int numcell = 0;
    int num_new_rgns = 0;  //number of tmpcalcID regions coming from oldrgn
    int ipttmp = pos_smrgnid;
    int jpt = j+old_neg_smrgnid;
    for (int i=(new_size_lookups-pos_smrgnid); i<new_size_lookups; i++) {  
      if (track_rid_matrix[i][j] > 0) num_new_rgns++;
      if (track_rid_matrix[i][j] > numcell) {
        numcell = track_rid_matrix[i][j]; 
        ipttmp = i+neg_smrgnid;
      }
    }
    pt_oldid_num_newrgns[jpt] = num_new_rgns;  //will contain number of new rgn from old rgn
    if (numcell > pt_smlookupvchk[ipttmp]) {
      pt_smlookupvchk[ipttmp] = numcell;
      pt_smlookup_oldid[ipttmp] = jpt;
    }
  }
  /* for negative regions */
  for (int j=0; j<(-old_neg_smrgnid); j++) {//j is browsing through old smrgnid with index
    int numcell = 0;
    int num_new_rgns = 0;  //number of tmpcalcID regions coming from oldrgn
    int ipttmp = neg_smrgnid;
    int jpt = j+old_neg_smrgnid;
    for (int i=0; i<(-neg_smrgnid); i++) {   //i is browsing through new smrgnids with index 
      if (track_rid_matrix[i][j] > 0) num_new_rgns++;
      if (track_rid_matrix[i][j] > numcell) {
        numcell = track_rid_matrix[i][j]; 
        ipttmp = i+neg_smrgnid;
      }
    }
    pt_oldid_num_newrgns[jpt] = num_new_rgns;  //will contain number of new rgn from old rgn
    if (numcell > pt_smlookupvchk[ipttmp]) {
      pt_smlookupvchk[ipttmp] = numcell;
      pt_smlookup_oldid[ipttmp] = jpt;
    }
  }
  boil::timer.stop("id_new->old");
 
  boil::timer.start("id_hiding");
  for (int j=0; j<old_size_lookups; j++) { //browse all old in track_rid_matrix
    int jid = j+old_neg_smrgnid;
    if (jid==0) continue;
    int totcell = 0;
    bool maybe_hiding = false;
    for (int i=0; i<new_size_lookups; i++) { //browse all new 
      int iid = i+neg_smrgnid;
      if (iid==0) continue;
      else if (track_rid_matrix[i][j] == 0) continue;
      else {  // if track_rid_matrix[i][j] > 0
        if (iid*jid > 0) { //same phase newID matches some of oldID
          maybe_hiding = false; //no hiding or removing
          break; 
        }
        else if (iid*jid < 0) {  //other phase newID matches some oldID 
          maybe_hiding = true;  //could be going into hiding
          totcell += track_rid_matrix[i][j];
        }
      }
    }
    if (maybe_hiding) { //other phase newID matches all oldID 
      Region & orgn = getregion(jid);
      if (totcell <= size_smallrgn) { //small and hiding !!
        if (!boil::cart.iam()) {
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<jid
                <<" cvol "<<totcell<<" xyz "<<orgn.x()<<" "<<orgn.y()
                <<" "<<orgn.z()<<" has gone into hiding"<<std::endl;
        }
        orgn.hiding(true);
//        need2plot = true;
      }
      else { //larger, entirely replaced with other phase
        if (!boil::cart.iam()) {
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<jid
                <<" cvol "<<totcell<<" xyz "<<orgn.x()<<" "<<orgn.y()
                <<" "<<orgn.z()<<" larger region erased"<<std::endl;
        }
        erase_region(jid);
      }
    }
  }
  boil::timer.stop("id_hiding");

  boil::timer.start("id_split");
  /* for all regions */
  for (int jo=old_neg_smrgnid; jo<=old_pos_smrgnid; jo++) {
    if (jo==0) continue;
    if (pt_oldid_num_newrgns[jo] > 1) {  //a split occured
      Region & orgn = getregion(jo);
      if (!boil::cart.iam()) {
#ifdef VERBOSE
        outrgn<<"Info3:t "<<time->current_time()<<" old ID "<<jo<<" cvol "<<orgn.cellvol()
              <<" xyz "<<orgn.x()<<" "<<orgn.y()<<" "<<orgn.z()
              <<" uvw "<<orgn.u()<<" "<<orgn.v()<<" "<<orgn.w()
              <<" splits into "<<pt_oldid_num_newrgns[jo]<<" regions: ";
#else
        outrgn<<"Info3:t "<<time->current_time()<<" old ID "<<jo
              <<" splits into "<<pt_oldid_num_newrgns[jo]<<" regions: ";
#endif
      }
      for (int in=neg_smrgnid; in<=pos_smrgnid; in++) {
        if (in==0) continue;
        if (track_rid_matrix[in-neg_smrgnid][jo-old_neg_smrgnid] > 0) { // tmpcalcID from split
          if (in*jo > 0) { //same phase newID matches some of oldID
            if (pt_smlookup_oldid[in] == 0) { //smaller id from split
              int newid = 0;
              if (in > 0) newid = get_new_pos_id(); //returns pos unsused id or new id
              else  newid = get_new_neg_id(); //returns neg unsused id or new id
              pt_smlookup_oldid[in] = newid; //give smaller split region new ID
              if (!boil::cart.iam()) {
#ifdef VERBOSE
                outrgn<<" new ID "<<newid<<" cvol "<<pt_smlookup_cvol[in]
                      <<" xyz "<<pt_flookup_x[in]<<" "<<pt_flookup_y[in]<<" "<<pt_flookup_z[in]
                      <<" uvw "<<pt_flookup_u[in]<<" "<<pt_flookup_v[in]<<" "<<pt_flookup_w[in];
#else
                outrgn<<" new ID "<<newid;
#endif
              }
            }
            else {  //larger id from split
              if (!boil::cart.iam()) {
#ifdef VERBOSE
                outrgn<<" cont ID "<<pt_smlookup_oldid[in]<<" cvol "<<pt_smlookup_cvol[in]
                      <<" xyz "<<pt_flookup_x[in]<<" "<<pt_flookup_y[in]<<" "<<pt_flookup_z[in]
                      <<" uvw "<<pt_flookup_u[in]<<" "<<pt_flookup_v[in]<<" "<<pt_flookup_w[in];
#else
                outrgn<<" cont ID "<<pt_smlookup_oldid[in];
#endif
              }
            }
          }
        }
      }
      outrgn<<std::endl;
    }
  }
  boil::timer.stop("id_split");

  boil::timer.start("id_output_hiding");
  for (int i=neg_smrgnid; i<=pos_smrgnid; i++) {
    if (i==0) continue;
    /* if 0 then tmpcalcID matches no oldID */
    if (pt_smlookup_oldid[i] == 0) { //for tmpcalcID i
      int cellvol = pt_smlookup_cvol[i];
      if (cellvol <= size_smallrgn) {  //if small vol
        int matchid = hidden_rgnid(pt_flookup_x[i], pt_flookup_y[i], pt_flookup_z[i]);
        if (matchid*i > 0) { //if true (nonzero) than hiding rgn with similar xyz, and must be same phase
          pt_smlookup_oldid[i] = matchid;
          if (!boil::cart.iam()) {
#ifdef VERBOSE
            outrgn<<"Info3:t "<<time->current_time()<<" ID "<<matchid<<" cvol "<<cellvol
                <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
                <<" uvw "<<pt_flookup_u[i]<<" "<<pt_flookup_v[i]<<" "<<pt_flookup_w[i]
                <<" is coming out of hiding"<<std::endl;
#else
            outrgn<<"Info3:t "<<time->current_time()<<" ID "<<matchid<<" cvol "<<cellvol
                <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
                <<" is coming out of hiding"<<std::endl;
#endif
          }
        }
        else { //small, not from split, no hiding rgn with similar xyz
          int newid = 0;
          if (i > 0) newid = get_new_pos_id(); //returns pos unsused id or new id
          else  newid = get_new_neg_id(); //returns neg unsused id or new id
          pt_smlookup_oldid[i] = newid;
          if (!boil::cart.iam()) {
#ifdef VERBOSE
            outrgn<<"Info3:t "<<time->current_time()<<" ID "<<newid<<" cvol "<<cellvol
                <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
                <<" uvw "<<pt_flookup_u[i]<<" "<<pt_flookup_v[i]<<" "<<pt_flookup_w[i]
                <<" is a small new ID"<<std::endl;
#else
            outrgn<<"Info3:t "<<time->current_time()<<" ID "<<newid<<" cvol "<<cellvol
                <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
                <<" is a small new ID"<<std::endl;
#endif
          }
        }
      } //end if small vol
      else { // if large  (vol > size_smallrgn)
        int newid = 0;
        if (i > 0) newid = get_new_pos_id(); //returns pos unsused id or new id
        else  newid = get_new_neg_id(); //returns neg unsused id or new id
        pt_smlookup_oldid[i] = newid;
        if (!boil::cart.iam()) {
#ifdef VERBOSE
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<newid<<" cvol "<<cellvol
              <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
              <<" uvw "<<pt_flookup_u[i]<<" "<<pt_flookup_v[i]<<" "<<pt_flookup_w[i]
              <<" is a large new ID"<<std::endl;
#else
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<newid<<" cvol "<<cellvol
              <<" xyz "<<pt_flookup_x[i]<<" "<<pt_flookup_y[i]<<" "<<pt_flookup_z[i]
              <<" is a large new ID"<<std::endl;
#endif
        }
      }
    }
  }
  boil::timer.stop("id_output_hiding");

  boil::timer.start("id_output_merging");
  /* for positive regions */
  for (int i=(new_size_lookups-pos_smrgnid); i<new_size_lookups; i++) {
    for (int j=(old_size_lookups-old_pos_smrgnid); j<old_size_lookups; j++) {
      if (track_rid_matrix[i][j] > 0) {
        int pti = i+neg_smrgnid; //tmpcalcID
        int ptj = j+old_neg_smrgnid;
        if (pt_smlookup_oldid[pti] != ptj) { //if tmpcalcID not using the oldID
          bool oldidused = false;
          for (int ii=neg_smrgnid; ii<=pos_smrgnid; ii++) {
            if (pt_smlookup_oldid[ii] == ptj) {oldidused=true; break;}
          }
          if (oldidused == false) { //no region using old id ptj
            Region & orgn = getregion(ptj); //returns old region not used
            int uid = pt_smlookup_oldid[pti];
            Region & urgn = getregion(uid);  //returns old region used by tmpcalcID
            if (!boil::cart.iam()) {
#ifdef VERBOSE
              outrgn<<"Info3:t "<<time->current_time()<<" ID "<<ptj<<" cvol "<<orgn.cellvol()
                    <<" xyz "<<orgn.x()<<" "<<orgn.y()<<" "<<orgn.z()
                    <<" uvw "<<orgn.u()<<" "<<orgn.v()<<" "<<orgn.w()
                    <<" and ID "<<uid<<" cvol "<<urgn.cellvol()
                    <<" xyz "<<urgn.x()<<" "<<urgn.y()<<" "<<urgn.z()
                    <<" uvw "<<urgn.u()<<" "<<urgn.v()<<" "<<urgn.w()
                    <<" merge into cont ID "<<uid<<" cvol "<<pt_smlookup_cvol[pti]
                    <<" xyz "<<pt_flookup_x[pti]<<" "<<pt_flookup_y[pti]<<" "<<pt_flookup_z[pti]
                    <<" uvw "<<pt_flookup_u[pti]<<" "<<pt_flookup_v[pti]<<" "<<pt_flookup_w[pti]
                    <<" ; "<<ptj<<" is now erased"<<std::endl;
#else
              outrgn<<"Info3:t "<<time->current_time()<<" ID "<<ptj
                    <<" and ID "<<uid
                    <<" merge into ID "<<uid
                    <<" ; "<<ptj<<" is erased"<<std::endl;
#endif
            }
            erase_region(ptj);
//            need2plot = true;
  } } } } }
  /* for negative regions */
  for (int i=0; i<(-neg_smrgnid); i++) {
    for (int j=0; j<(-old_neg_smrgnid); j++) {
      if (track_rid_matrix[i][j] > 0) {
        int pti = i+neg_smrgnid; //convert back to IDs for pointers
        int ptj = j+old_neg_smrgnid;
        if (pt_smlookup_oldid[pti] != ptj) { //if tmpcalcID not the oldID
          bool oldidused = false;
          for (int ii=neg_smrgnid; ii<=pos_smrgnid; ii++) {
            if (pt_smlookup_oldid[ii] == ptj){oldidused=true; break;}
          }
          if (oldidused == false) { //no region using old id ptj
            Region & orgn = getregion(ptj); //returns old region not used
            int uid = pt_smlookup_oldid[pti];
            Region & urgn = getregion(uid);  //returns old region used by tmpcalcID
            if (!boil::cart.iam()) {
#ifdef VERBOSE
              outrgn<<"Info3:t "<<time->current_time()<<" ID "<<ptj<<" cvol "<<orgn.cellvol()
                    <<" xyz "<<orgn.x()<<" "<<orgn.y()<<" "<<orgn.z()
                    <<" uvw "<<orgn.u()<<" "<<orgn.v()<<" "<<orgn.w()
                    <<" and ID "<<uid<<" cvol "<<urgn.cellvol()
                    <<" xyz "<<urgn.x()<<" "<<urgn.y()<<" "<<urgn.z()
                    <<" uvw "<<urgn.u()<<" "<<urgn.v()<<" "<<urgn.w()
                    <<" merge into cont ID "<<uid<<" cvol "<<pt_smlookup_cvol[pti]
                    <<" xyz "<<pt_flookup_x[pti]<<" "<<pt_flookup_y[pti]<<" "<<pt_flookup_z[pti]
                    <<" uvw "<<pt_flookup_u[pti]<<" "<<pt_flookup_v[pti]<<" "<<pt_flookup_w[pti]
                    <<" ; "<<ptj<<" is now erased"<<std::endl;
#else
              outrgn<<"Info3:t "<<time->current_time()<<" ID "<<ptj
                    <<" and ID "<<uid
                    <<" merge into ID "<<uid
                    <<" ; "<<ptj<<" is erased"<<std::endl;
#endif
            }
            erase_region(ptj);
//            need2plot = true;
  } } } } }

  boil::timer.stop("id_output_merging");

  for (int i=neg_smrgnid; i<=pos_smrgnid; i++) { //for all regions
    if (i==0) continue;  //do not create a new region with id=0
    int oldid = pt_smlookup_oldid[i];
    Region & updatergn = getregion(oldid);
    updatergn.cellvol(pt_smlookup_cvol[i]);
    updatergn.xyz(pt_flookup_x[i], pt_flookup_y[i], pt_flookup_z[i]); //set new pos
    updatergn.uvw(pt_flookup_u[i], pt_flookup_v[i], pt_flookup_w[i]); //set new vel
  }
  
  /* erase region that is not hiding and not used */
  for (int i=0; i<m_vectrgns.size(); i++) {
    if (m_vectrgns[i].hiding() == false) {
      int sid = m_vectrgns[i].id();
      bool rgnused = false;
      for (int itmp=neg_smrgnid; itmp<=pos_smrgnid; itmp++) { //browse through used IDs
        if (itmp==0) continue;
        else if (sid==pt_smlookup_oldid[itmp]) {
          rgnused = true;
          break;
        }
      } 
      if (rgnused == false) {
        Region & trgn = m_vectrgns[i];
        if (!boil::cart.iam()) {
#ifdef VERBOSE
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<sid<<" cvol "<<trgn.cellvol()
                <<" xyz "<<trgn.x()<<" "<<trgn.y()<<" "<<trgn.z()
                <<" not used so erased"<<std::endl;
#else
          outrgn<<"Info3:t "<<time->current_time()<<" ID "<<sid
                <<" not used so erased"<<std::endl;
#endif
        }
        erase_region(sid);
//        need2plot = true;
      }
    }
  }
        
  /* erase hidden regions that have been in hiding too long */
  mng_hidden_rgns();

  /* update rgnid field with old regions id and new added region ids  */
  for_vijk(rgnid,i,j,k) {
    int rid = int(rgnid[i][j][k]);
    rgnid[i][j][k] = pt_smlookup_oldid[rid];
  }

  /* output region vol, com (xyz), <u> (uvw) every out_rgn_info_freq */
  if (time->current_step() % out_rgn_info_freq == 0)
    out_region_info();

  /* never plot */
//  need2plot = false;  //to quickly toggle need2plot

  /*plot dbg */
  if (need2plot) {
    if (!boil::cart.iam()) {
      outrgn<<"Info3:t "<<time->current_time()<<" Plotting: c-rgn-rgnold "<<pltinc<<std::endl;
    }
    boil::plot->plot(c, rgnid, rgnid_old, "c-rgn-rgnold", pltinc);
    pltinc++;
  }
  /*plot dbg */

  boil::timer.stop("id_rgns_mn");

}
