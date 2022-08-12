/***************************************************************************//**
 * cangle
*******************************************************************************/
inline real cangle(const int i, const int j, const int k) const {
  if(!cangle_variable)
    return cangle0;
  else
    return cangle_func(i,j,k)/180.0*boil::pi;
}

/* getter for cangle */
inline real get_cangle() const { return(cangle0/boil::pi*180.0);};

/* setter for cangle */
inline void set_cangle(const real r) {
  cangle_variable = false;
  cangle0=r/180.0*boil::pi;
  boil::oout<<"VOF::set_cangle: cangle= "<<r<<"\n";
}

inline void set_cangle(const boil::func_ijk_real & f) {
  cangle_variable = true;
  cangle_func = f;
  boil::oout<<"VOF::cangle set variable"<<"\n";
}

/***************************************************************************//**
 * update at walls
*******************************************************************************/
inline void set_update_at_walls_variable(const boil::func_scijk_real & f) {
  update_at_walls_variable = true;
  update_at_walls_func = f;
  boil::oout<<"VOF::update_at_walls set variable"<<"\n";
}

inline void set_update_at_walls_variable(const bool st) {
  update_at_walls_variable = st;
  assert(!st); /* this should be used for turning variability off only */
  boil::oout<<"VOF::update_at_walls set variable "<<st<<"\n";
}

/***************************************************************************//**
 * wall curvature height construction
*******************************************************************************/
/* getter for dir */
inline Comp get_wall_curv_dir() const { return wall_curv_dir;};

/* setter for dir */
inline void set_wall_curv_dir(const Comp & m) {
  wall_curv_dir = m;
  if(!boil::cart.iam()) {
    boil::oout<<"VOF::wall_curv_dir: "<<m<<"\n";
  }
}

/***************************************************************************//**
 *  other
*******************************************************************************/
/* getter and setter for wall value tolerance */
inline real get_tol_wall() const { return tol_wall; }
inline void set_tol_wall(real tolnew) {
  tol_wall = tolnew;
  boil::oout<<"VOF: New wall value tolerance: "<<tol_wall<<boil::endl;
  return;
}

/* getter and setter for bulk_curv_method */
inline CurvMethod get_curv_method() const {return bulk_curv_method;}
void set_curv_method(const CurvMethod cm) {
  bulk_curv_method=cm;
  if(cm==CurvMethod::HF()) {
    boil::oout<<"VOF: height function is used for curvature calculation.\n";
  } else {
    boil::oout<<"Curv method should be HF.\n";
    boil::oout<<"Exiting."<<boil::endl;
    exit(0);
  }
}

/* getter and setter for advect_method */
inline AdvectionMethod get_advection_method() const {return advect_method;}
void set_advection_method(const AdvectionMethod am) {
  advect_method=am;
  if(!boil::cart.iam())
    boil::oout<<"VOF:::advection method: "<<am<<boil::endl;
}

/* setter for limit_color */
inline void set_limit_color(const bool b) {
  limit_color=b;
  boil::oout<<"set_limit_color= "<<b<<"\n";
}
/* getter for limit_color */
inline bool get_limit_color() const { return(limit_color);};

/* setter for subgrid method */
inline void set_subgrid_method(const SubgridMethod b) {
  subgrid_method=b;
  if(!boil::cart.iam())
    boil::oout<<"Subgrid method: "<<subgrid_method<<boil::endl;
}
/* getter for subgrid method */
inline SubgridMethod get_subgrid_method() const { return(subgrid_method);};

/* setter for use_interp */
inline void set_use_interp(const bool b) {
  use_interp=b;
  boil::oout<<"set_use_interp= "<<b<<"\n";
}
/* getter for use_interp */
inline bool get_use_interp() const { return(use_interp);};

/* min and max of color function in fluid domain */
inline real minval() const {return minclr;}
inline real maxval() const {return maxclr;}
inline void set_minval(real r) {minclr=r;}
inline void set_maxval(real r) {maxclr=r;}

/* setter for normal vector method */
void set_normal_vector_method_advance(const NormMethod nm) {
  norm_method_advance = nm;
  if(!boil::cart.iam())
    boil::oout<<"Normal vector method for advance: "<<nm<<boil::endl;
  if(nm==NormMethod::ElviraYZ()) {
    mcomp_for_elvira = Comp::i();
  } else if(nm==NormMethod::ElviraXZ()) {
    mcomp_for_elvira = Comp::j();
  } else if(nm==NormMethod::ElviraXY()) {
    mcomp_for_elvira = Comp::k();
  } else {
    mcomp_for_elvira = Comp::undefined();
  }
}
void set_normal_vector_method_curvature(const NormMethod nm) {
  norm_method_curvature = nm;
  if(!boil::cart.iam())
    boil::oout<<"Normal vector method for curvature: "<<nm<<boil::endl;
  if(nm==NormMethod::ElviraYZ()) {
    mcomp_for_elvira = Comp::i();
  } else if(nm==NormMethod::ElviraXZ()) {
    mcomp_for_elvira = Comp::j();
  } else if(nm==NormMethod::ElviraXY()) {
    mcomp_for_elvira = Comp::k();
  } else {
    mcomp_for_elvira = Comp::undefined();
  }
}
void set_normal_vector_method_all(const NormMethod nm) {
  norm_method_advance = nm;
  norm_method_curvature = nm;
  if(!boil::cart.iam())
    boil::oout<<"Normal vector method: "<<nm<<boil::endl;
  if(nm==NormMethod::ElviraYZ()) {
    mcomp_for_elvira = Comp::i();
  } else if(nm==NormMethod::ElviraXZ()) {
    mcomp_for_elvira = Comp::j();
  } else if(nm==NormMethod::ElviraXY()) {
    mcomp_for_elvira = Comp::k();
  } else {
    mcomp_for_elvira = Comp::undefined();
  }
}

/* getter for normal vector method */
inline NormMethod get_normal_vector_method_advance() const {
  return norm_method_advance;
}
inline NormMethod get_normal_vector_method_curvature() const {
  return norm_method_curvature;
}

#if 1
/* setter for near-wall curvature method */
void set_wall_curv_method(const CurvMethod wcm) {
  wall_curv_method = wcm;
  boil::oout<<"vof:set wall_curv_method= "<<wall_curv_method<<"\n";
}
CurvMethod get_wall_curv_method() {
  return wall_curv_method;
}
#endif

#if 0
/* setter for near-wall curvature method */
void set_wall_curv_method(const CurvMethod wcm,
                          const Sign sig = Sign::undefined(),
                          const real cgl = -1.,
                          const int nfilm_crit = -1) {
  wall_curv_method = wcm;
  if(!boil::cart.iam())
    boil::oout<<"Wall curvature method: "<<wcm<<boil::endl;
  if(wcm==CurvMethod::HFparallelXZ()) {//||wcm==CurvMethod::HFmixedXZ()) {
    if       (sig==Sign::pos()) {
      mult_wall =  1;
    } else if(sig==Sign::neg()) {
      mult_wall = -1;
    } else {
      boil::oout<<"Phase at the origin must be specified! "
                <<"Sign-neg: phi(orig)<0.5 and vice versa. "
                <<"Exiting."
                <<boil::endl;
      exit(0);
    }
  } else {
    if(sig!=Sign::undefined()) {
      if(sig==Sign::pos()) {
        mult_wall =  1;
      } else {
        mult_wall = -1;
      }
    }
  }
#if 0
  if(wcm==CurvMethod::HFmixedXZ()) {
    if       (nfilm_crit>0) {
      Nfilm_crit = nfilm_crit;
    } else if(cgl>0.) {
      /* estimate */
      real cangrad = cgl*boil::pi/180.;
      Nfilm_crit = std::max(4,abs(int(cos(cangrad)/sin(cangrad))));
      boil::oout<<"Wall critical film length: "<<Nfilm_crit<<boil::endl;
    } else {
      Nfilm_crit = boil::unint;
    }
  }
#endif

  if(wcm==CurvMethod::HFparallelXZ()) {
    if(cgl>0.) {
      detachment_model.set_detachment_params(cgl);
    }
  }
  return;
}

/* getter for near-wall curvature method */
inline CurvMethod get_wall_curv_method() const { return wall_curv_method;};

/* getter for critical film length */
inline int get_critical_film_length() const { return Nfilm_crit;};
#endif

/* setter for pressure extrapolation */
void set_pressure_extrapolation_parameters(const bool eflag, const int eiter) {
  store_pressure_extrap = eflag;
  niter_pressure_extrap = eiter;
  boil::oout<<"set_pressure_extrapolation_parameters: store_pressure= "<<eflag<<" ; "
                <<"maxiter= "<<eiter<<"\n";
  return;
}

/* getter for pressure extrapolation */
inline bool get_store_pressure_extrapolation_flag() const { return store_pressure_extrap;}
inline int get_pressure_extrapolation_maxiter() const { return niter_pressure_extrap;}

/* bounded color */
inline real bounded_color(const real cval) const {
  return std::max(0.0,std::min(1.0,cval));
}

/* setter for extrapolate_ib */
void set_extrapolate_ib(const bool b) {
  extrapolate_ib=b;
  boil::oout<<"cipcsl2:extrapolate_ib= "<<b<<"\n";
  boil::oout<<"        true:  extrapolate color function into solid \n";
  boil::oout<<"        false: no extrapolate color function into solid \n";
}
bool get_extrapolate_ib() { return (extrapolate_ib); };

