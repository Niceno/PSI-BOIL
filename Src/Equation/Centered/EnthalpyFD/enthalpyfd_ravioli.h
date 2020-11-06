////////////////
//            //
//  Stencil   //
//            //
////////////////
/* this is a ravioli struct for stencil  */
struct StencilPoint {
  public:
    StencilPoint(int i = 0, real v = 0., real p = 0.) :
      idx(i), val(v), pos(p) {};

    int idx;
    real val, pos;
};
