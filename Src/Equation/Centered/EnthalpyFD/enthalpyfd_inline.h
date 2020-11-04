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
