    inline int get_gradt_accuracy_order() const { return ao_conv.eval(); }
    inline void set_gradt_accuracy_order(const int ao) {
      ao_conv = AccuracyOrder(ao);
      if(!boil::cart.iam())
        boil::oout<<"EnthalpyFD::gradt_accuracy_order= "
                  <<ao<<"\n";
    }
    inline void set_gradt_accuracy_order(const AccuracyOrder ao) {
      ao_conv = ao;
      if(!boil::cart.iam())
        boil::oout<<"EnthalpyFD::gradt_accuracy_order= "
                  <<ao<<"\n";
    }

    inline bool get_no_solid_acceleration() const
      { return accelerated_no_solid; }
    inline void set_no_solid_acceleration(const bool flag) {
      accelerated_no_solid = flag;
      boil::oout<<"EnthalpyFD::no_solid_acceleration= "
                <<accelerated_no_solid<<"\n";
    }

    inline bool stencil_min(const Comp & m, 
                            const int i, const int j, const int k) const {
      return (m==Comp::i())*(i==si() && bflag_struct.iminc)
            +(m==Comp::j())*(j==sj() && bflag_struct.jminc)
            +(m==Comp::k())*(k==sk() && bflag_struct.kminc);
    }

    inline bool stencil_max(const Comp & m, 
                            const int i, const int j, const int k) const {
      return (m==Comp::i())*(i==ei() && bflag_struct.imaxc)
            +(m==Comp::j())*(j==ej() && bflag_struct.jmaxc)
            +(m==Comp::k())*(k==ek() && bflag_struct.kmaxc);
    }

    inline const Vector & flx_liq() const { return flux_liq; }
    inline const Vector & flx_gas() const { return flux_gas; }
