    inline int get_extrapolation_iters() const { return mmax_ext; }
    inline real get_extrapolation_tol() const  { return tol_ext; }
    inline void set_extrapolation_params(const int mnew, const real tolnew) {
      mmax_ext = mnew;
      tol_ext = tolnew;
      boil::oout<<"Topology::extrapolationparams: "<<mnew<<" "<<tolnew<<"\n";
    }

    inline real get_close_to_cc() const { return close_to_cc; }
    inline void set_close_to_cc(const real ccc_new) {
      close_to_cc = ccc_new;
      boil::oout<<"Topology::close_to_cc: "<<close_to_cc<<"\n";
    }

    inline bool above_interface(const real c)
      const { return c>=clrsurf; }
    inline bool below_interface(const real c)
      const { return c<clrsurf; }
    inline bool above_interface(const int i, const int j, const int k)
      const { return (*clr)[i][j][k]>=clrsurf; }
    inline bool above_interface_old(const int i, const int j, const int k)
      const { return clrold[i][j][k]>=clrsurf; }
    inline bool below_interface(const int i, const int j, const int k)
      const { return (*clr)[i][j][k]<clrsurf; }
    inline bool below_interface_old(const int i, const int j, const int k)
      const { return clrold[i][j][k]<clrsurf; }

    inline bool above_interface(const int i, const int j, const int k,
                                const Old old) const {
      return (old==Old::yes) ?
             above_interface_old(i,j,k) :
             above_interface(i,j,k);
    }
 
    inline bool below_interface(const int i, const int j, const int k,
                                const Old old) const {
      return (old==Old::yes) ?
             below_interface_old(i,j,k) :
             below_interface(i,j,k);
    }

    /* <0: below interface, >0 above interface */
    inline Sign sign_interface(const real c) const {
      return above_interface(c) ? Sign::pos() : Sign::neg();
    }
    inline Sign sign_interface(const int i, const int j, const int k,
                               const Old old) const { 
      return (old==Old::yes) ?
             sign_interface_old(i,j,k) : sign_interface(i,j,k);
    }
    inline Sign sign_interface(const int i, const int j, const int k) 
      const { return sign_interface((*clr)[i][j][k]); }
    inline Sign sign_interface_old(const int i, const int j, const int k) 
      const { return sign_interface(clrold[i][j][k]); }

/***************************************************************************//*
 *  distances
******************************************************************************/
inline real distance_int(const Sign dir, const Comp & m,
                         const int i, const int j, const int k,
                         Sign & cell_marker, const Old old) const {
  return (old==Old::yes) ?
         distance_int_old(dir,m,i,j,k,cell_marker) :
         distance_int(dir,m,i,j,k,cell_marker);
}

inline real distance_int_x(const Sign dir,
                           const int i, const int j, const int k,
                           Sign & cell_marker, const Old old) const {
  return (old==Old::yes) ?
         distance_int_x_old(dir,i,j,k,cell_marker) :
         distance_int_x(dir,i,j,k,cell_marker);

}
inline real distance_int_y(const Sign dir,
                           const int i, const int j, const int k,
                           Sign & cell_marker, const Old old) const {
  return (old==Old::yes) ?
         distance_int_y_old(dir,i,j,k,cell_marker) :
         distance_int_y(dir,i,j,k,cell_marker);

}
inline real distance_int_z(const Sign dir,
                           const int i, const int j, const int k,
                           Sign & cell_marker, const Old old) const {
  return (old==Old::yes) ?
         distance_int_z_old(dir,i,j,k,cell_marker) :
         distance_int_z(dir,i,j,k,cell_marker);
}

/***************************************************************************//*
 *  differences; order = number of points used (not acc order)
******************************************************************************/
    /* zeroth difference = interpolation */
    real nth_order_nth(const std::vector<StencilPoint> & stencil,
                       const AccuracyOrder & order_der,
                       const AccuracyOrder & order_dif) const {
      switch(order_der.eval()) {
        case 0 :
          return nth_order_zeroth(stencil,order_dif);
        case 1 :
          return nth_order_first(stencil,order_dif);
        case 2 :
          return nth_order_second(stencil,order_dif);
        case 3 :
          return nth_order_third(stencil,order_dif);
        default :
          boil::aout<<"Topology: unrecognised difference requested. Exiting."
                    <<boil::endl;
          exit(0);
      }
      return 0.0;
    }

    /* getter for front_minmax */
    inline real get_xminft() const { return(xminft);};
    inline real get_xmaxft() const { return(xmaxft);};
    inline real get_yminft() const { return(yminft);};
    inline real get_ymaxft() const { return(ymaxft);};
    inline real get_zminft() const { return(zminft);};
    inline real get_zmaxft() const { return(zmaxft);};

    /* getter for area */
    inline real get_totarea() const { 
      real are(0.);
      for_vijk(get_adens(),i,j,k)
        are += get_adens()[i][j][k]*get_adens().dV(i,j,k);
      boil::cart.sum_real(&are);
      return are;
    }
    inline real get_area(const int i, const int j, const int k) const {
      return get_adens()[i][j][k]*get_adens().dV(i,j,k);
    }

    /* references */
    Scalar & get_vf()    { return *vf; }
    Scalar & get_clr()   { return *clr; }
    Scalar & get_nx()    { return *nx; }
    Scalar & get_ny()    { return *ny; }
    Scalar & get_nz()    { return *nz; }
    Scalar & get_adens() { return *adens; }
    Vector & get_fs()    { return *fs; }
    ScalarInt & get_iflag() { return *iflag; }

    const Scalar & get_vf() const    { return *vf; }
    const Scalar & get_clr() const   { return *clr; }
    const Scalar & get_nx() const    { return *nx; }
    const Scalar & get_ny() const    { return *ny; }
    const Scalar & get_nz() const    { return *nz; }
    const Scalar & get_adens() const { return *adens; }
    const Vector & get_fs() const    { return *fs; }
    const ScalarInt & get_iflag() const { return *iflag; }

    const Domain * domain() const {return clrold.domain();}
