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
                    const Scalar * diff_eddy) const {
  return diff_eddy ? val_lambdal + (*diff_eddy)[i][j][k]*val_cpl/val_rhol/turbP
                   : val_lambdal;
}
inline real lambdav(const int i, const int j, const int k,
                    const Scalar * diff_eddy) const {
  return diff_eddy ? val_lambdav + (*diff_eddy)[i][j][k]*val_cpv/val_rhov/turbP
                   : val_lambdav;
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
                      const int i, const int j, const int k) const{

  return topo->interface(dir,m,i,j,k);
}

inline bool interface_old(const Sign dir, const Comp m,
                          const int i, const int j, const int k) const {

  return topo->interface_old(dir,m,i,j,k);
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
