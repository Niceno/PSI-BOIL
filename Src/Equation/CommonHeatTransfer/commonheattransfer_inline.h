/***************************************************************************//**
 *  Material properties 
******************************************************************************/
inline real cpl(const int i, const int j, const int k) const {
  return val_cpl;
}
inline real cpv(const int i, const int j, const int k) const {
  return val_cpv;
}
inline real rhol(const int i, const int j, const int k) const {
  return val_rhol;
}
inline real rhov(const int i, const int j, const int k) const {
  return val_rhov;
}
inline real lambdal(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const {
  return diff_eddy ? val_lambdal + (*diff_eddy)[i][j][k]*val_cpl/val_rhol/turbP
                   : val_lambdal;
}
inline real lambdav(const int i, const int j, const int k,
                    const Scalar * diff_eddy = NULL) const {
  return diff_eddy ? val_lambdav + (*diff_eddy)[i][j][k]*val_cpv/val_rhov/turbP
                   : val_lambdav;
}

/***************************************************************************//**
 * Heat transfer resistance and wall source
******************************************************************************/
inline real int_resistance_liq(const int i, const int j, const int k) const {
  if(!int_resistance_liq_variable)
    return int_resistance_liq_val;
  else
    return int_resistance_func(i,j,k);
}

#if 0
inline real int_resistance_vap(const int i, const int j, const int k) const {
  return int_resistance_vap_val;
}
#endif

inline real wall_resistance(const int i, const int j, const int k) const {
  if(!wall_resistance_variable)
    return wall_resistance_val;
  else
    return wall_resistance_func(i,j,k);
}

inline real dirac_wall_source(const int i, const int j, const int k) const {
  if(!dirac_wall_source_variable)
    return dirac_wall_source_val;
  else
    return dirac_wall_source_func(i,j,k);
}

/***************************************************************************//**
*  \brief calculate temperature at node point
*             len_s         len_f
*             lam_s         lam_f
*         *-------------*------------*
*       tmp_s       tmp_node        tmp_f
*******************************************************************************/
inline real temperature_node(
             const real len_s, const real lam_s, 
             const real tmp_s, const real len_f,
             const real lam_f, const real tmp_f) const {
  //return (len_f*lam_s*tmp_s + len_s*lam_f*tmp_f)/(len_f*lam_s + len_s*lam_f);
  return temperature_node(0.,len_s/lam_s, tmp_s, len_f/lam_f, tmp_f);
}

inline real temperature_node(
             const real R_s, const real tmp_s,
             const real R_f, const real tmp_f) const {
  //return (R_f*tmp_s+R_s*tmp_f)/(R_f+R_s);
  return temperature_node(0.,R_s,tmp_s,R_f,tmp_f);
}

/* with a source */
inline real temperature_node(const real Q,
             const real R_s, const real tmp_s,
             const real R_f, const real tmp_f) const {
  return (R_f*R_s*Q+R_f*tmp_s+R_s*tmp_f)/(R_f+R_s);
}

/***************************************************************************//**
 *  Checks if the given cell is at (or next to) an interface
******************************************************************************/
inline bool interface(const Sign dir, const Comp m,
                      const int i, const int j, const int k, const Old old)
                      const {

  return (old==Old::yes) ? 
         topo->interface_old(dir,m,i,j,k) : topo->interface(dir,m,i,j,k);
}

inline bool interface(const int i, const int j, const int k) const {

  return topo->interface(i,j,k);
}

inline bool above_interface(const int i, const int j, const int k,
                            const Old old) const {
  return (old==Old::yes) ?
         topo->above_interface_old(i,j,k) : topo->above_interface(i,j,k);
}

inline bool below_interface(const int i, const int j, const int k,
                            const Old old) const {
  return (old==Old::yes) ?
         topo->below_interface_old(i,j,k) : topo->below_interface(i,j,k);
}

/***************************************************************************//**
 *  Call to tifmodel
******************************************************************************/
inline real Tint(const int dir, const Comp &mcomp, real frac,
                 const int i, const int j, const int k) const {
  return tifmodel.Tint(dir,mcomp,frac,i,j,k);
}

inline real Tint_old(const int dir, const Comp &mcomp, real frac,
                     const int i, const int j, const int k) const {
  return tifmodel.Tint_old(dir,mcomp,frac,i,j,k);
}

inline real Tint(const int i, const int j, const int k, const Old old) const {
  return (old==Old::yes) ?
         tifmodel.Tint_old(i,j,k) : tifmodel.Tint(i,j,k);
}

inline real Tint(const int i, const int j, const int k) const {
  return tifmodel.Tint(i,j,k);
}

inline real Tint_old(const int i, const int j, const int k) const {
  return tifmodel.Tint_old(i,j,k);
}

/***************************************************************************//**
 *  other
******************************************************************************/
inline real get_turbP() const { return turbP; }
inline void set_turbP(const real a) {
  turbP = a;
  boil::oout<<"CommonHeatTransfer::turbP= "<<turbP<<"\n";
}

/* units W/m2 */
inline real get_dirac_wall_source() const { return dirac_wall_source_val; }
inline void set_dirac_wall_source(const real a) {
  dirac_wall_source_variable = false;
  dirac_wall_source_val = a;
  boil::oout<<"CommonHeatTransfer::dirac_wall_source= "
            <<dirac_wall_source_val<<"\n";
}
inline void set_dirac_wall_source(const boil::func_ijk_real & f) {
  dirac_wall_source_variable = true;
  dirac_wall_source_func = f;
  boil::oout<<"CommonHeatTransfer::dirac_wall_source set variable"
            <<"\n";
}


/* units [K/(W/m^2)] = [m/(W/mK)] */
inline void set_int_resistance(const real re) {
  //int_resistance_vap_val = re;
  int_resistance_liq_val = re;
  use_int_resist = true;
  boil::oout<<"CommonHeatTransfer::int_resistance= "<<re<<"\n";
}
inline bool use_int_resistance() const {
  return use_int_resist;
}

#if 0
inline real get_int_resistance_vap() const { 
  return int_resistance_vap_val;
}
inline void set_int_resistance_vap(const real re) {
  int_resistance_vap_val = re;
  boil::oout<<"CommonHeatTransfer::int_resistance_vap= "
            <<re<<"\n";
}
#endif

inline real get_int_resistance_liq() const { 
  return int_resistance_liq_val;
}
inline void set_int_resistance_liq(const real re) {
  int_resistance_liq_val = re;
  use_int_resist = true;
  int_resistance_liq_variable = false;
  boil::oout<<"CommonHeatTransfer::int_resistance_liq= "
            <<re<<"\n";
}
inline void set_int_resistance_liq(const boil::func_ijk_real & f) {
  int_resistance_liq_variable = true;
  int_resistance_func = f;
  boil::oout<<"CommonHeatTransfer::int_resistance_liq set variable"
            <<"\n";
}

inline real get_wall_resistance() const {
  return wall_resistance_val;
}
inline void set_wall_resistance(const real re) {
  wall_resistance_variable = false;
  wall_resistance_val = re;
  boil::oout<<"CommonHeatTransfer::wall_resistance= "<<re<<"\n";
}
inline void set_wall_resistance(const boil::func_ijk_real & f) {
  wall_resistance_variable = true;
  wall_resistance_func = f;
  boil::oout<<"CommonHeatTransfer::wall_resistance set variable"
            <<"\n";
}

/**************************/

inline void init(const Scalar * diff_eddy = NULL) { new_time_step(diff_eddy); }
inline void new_time_step(const Scalar * diff_eddy = NULL) { 
  if(solid())
    calculate_node_temperature(diff_eddy); 
  return;
}
