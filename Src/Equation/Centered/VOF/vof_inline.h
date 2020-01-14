/* getter for front_minmax */
inline real get_xminft() { return(xminft);};
inline real get_xmaxft() { return(xmaxft);};
inline real get_yminft() { return(yminft);};
inline real get_ymaxft() { return(ymaxft);};
inline real get_zminft() { return(zminft);};
inline real get_zmaxft() { return(zmaxft);};

/* getter and setter for wall value tolerance */
inline real get_tol_wall() { return tol_wall; }
inline void set_tol_wall(real tolnew) {
  tol_wall = tolnew;
  boil::oout<<"VOF: New wall value tolerance: "<<tol_wall<<boil::endl;
  return;
}

/* getter and setter for bulk_curv_method */
inline CurvMethod get_curv_method() {return bulk_curv_method;}
void set_curv_method(const CurvMethod cm) {
  bulk_curv_method=cm;
  if(cm==CurvMethod::HF()){
    boil::oout<<"VOF: height function is used for curvature calculation.\n";
  } else if(cm==CurvMethod::DivNorm()){
    boil::oout<<"VOF: smoothed VOF is used for curvature calculation.\n";
  } else {
    boil::oout<<"Curv method should be HF or DivNorm.\n";
    boil::oout<<"Exiting."<<boil::endl;
    exit(0);
  }
}

/* setter for cangle */
inline void set_cangle(const real r) {
  cangle=r/180.0*boil::pi;
  boil::oout<<"set_cangle: cangle= "<<r<<"\n";
}
/* getter for cangle */
inline real get_cangle() { return(cangle/boil::pi*180.0);};

/* setter for limit_color */
inline void set_limit_color(const bool b) {
  limit_color=b;
  boil::oout<<"set_limit_color= "<<b<<"\n";
}
/* getter for limit_color */
inline bool get_limit_color() { return(limit_color);};

/* setter for subgrid method */
inline void set_subgrid_method(const SubgridMethod b) {
  subgrid_method=b;
  boil::oout<<"Subgrid method: "<<subgrid_method<<boil::endl;
}
/* getter for subgrid method */
inline SubgridMethod get_subgrid_method() { return(subgrid_method);};

/* setter for use_interp */
inline void set_use_interp(const bool b) {
  use_interp=b;
  boil::oout<<"set_use_interp= "<<b<<"\n";
}
/* getter for use_interp */
inline bool get_use_interp() { return(use_interp);};

/* min and max of color function in fluid domain */
inline real minval() {return minclr;}
inline real maxval() {return maxclr;}
void color_minmax();
inline void set_minval(real r) {minclr=r;}
inline void set_maxval(real r) {maxclr=r;}

/* setter for normal vector method */
void set_normal_vector_method_advance(const NormMethod nm) {
  norm_method_advance = nm;
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
inline NormMethod get_normal_vector_method_advance() {
  return norm_method_advance;
}
inline NormMethod get_normal_vector_method_curvature() {
  return norm_method_curvature;
}

/* setter for near-wall curvature method */
void set_wall_curv_method(const CurvMethod wcm,
                          const Sign sig = Sign::undefined(),
                          const real cangle = -1.) {
  wall_curv_method = wcm;
  boil::oout<<"Wall curvature method: "<<wcm<<boil::endl;
  if(wcm==CurvMethod::HFmixedXZ()) {
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

    if(cangle>0.) {
      detachment_model.set_detachment_params(cangle);
    }
  }
  return;
}

/* getter for near-wall curvature method */
inline CurvMethod get_wall_curv_method() { return wall_curv_method;};

/* setter for topoogy method */
void set_topo_method(const TopoMethod tpm) {
  topo_method = tpm;
  boil::oout<<"Topology method: "<<tpm<<boil::endl;
  return;
}

/* getter for topology method */
inline TopoMethod get_topo_method() { return topo_method;};

/* setter for pressure extrapolation */
void set_pressure_extrapolation_parameters(const bool eflag, const int eiter) {
  store_pressure_extrap = eflag;
  niter_pressure_extrap = eiter;
  boil::oout<<"set_pressure_extrapolation_parameters: store_pressure= "<<eflag<<" ; "
                <<"maxiter= "<<eiter<<"\n";
  return;
}

/* getter for pressure extrapolation */
inline bool get_store_pressure_extrapolation_flag() { return store_pressure_extrap;}
inline int get_pressure_extrapolation_maxiter() { return niter_pressure_extrap;}
