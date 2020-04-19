    inline real get_turbP() const { return turbP; }
    inline void set_turbP(const real a) {
      turbP = a;
      boil::oout<<"PhaseChange4::turbP= "<<turbP<<"\n";
    }
 
    inline bool get_accuracy_flag() const { return use_second_order_accuracy; }
    inline void set_second_order_accuracy(const bool flag) {
      use_second_order_accuracy = flag;
      boil::oout<<"PhaseChange4::use_second_order_accuracy= "
                <<use_second_order_accuracy<<"\n";
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
