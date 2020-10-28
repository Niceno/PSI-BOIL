    inline int get_accuracy_order() const { return accuracy_order.eval(); }
    inline void set_accuracy_order(const int ao) {
      accuracy_order = AccuracyOrder(ao);
      boil::oout<<"PhaseChange4::accuracy_order= "
                <<ao<<"\n";
    }
    inline void set_second_order_accuracy(const bool flag) {
      if(flag)
        accuracy_order = AccuracyOrder::Second();
      else
        accuracy_order = AccuracyOrder::First();
      boil::oout<<"PhaseChange4::accuracy_order= "
                <<accuracy_order<<"\n";
    }
 
    inline bool get_extrapolation_flag() const {
      return use_unconditional_extrapolation;
    }
    inline void set_unconditional_extrapolation(const bool flag) {
      use_unconditional_extrapolation = flag;
      boil::oout<<"PhaseChange4::use_unconditional_extrapolation= "
                <<use_unconditional_extrapolation<<"\n";
    }
    
    inline bool get_discard_flag() const {
      return discard_points_near_interface;
    }
    inline void set_discard_points_near_interface(const bool flag) {
      discard_points_near_interface = flag;
      boil::oout<<"PhaseChange4::discard_points_near_interface= "
                <<discard_points_near_interface<<"\n";
    }

    inline real get_smdot_pos() const {return smdot_pos;}
    inline real get_smdot_neg() const {return smdot_neg;}
    inline real get_smdot() const {return smdot_pos+smdot_neg;}
