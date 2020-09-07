#ifndef HTWM_H
#define HTWM_H

///////////////////////////////
//                           //
// Heat Transfer Wall Model  //
//                           //
///////////////////////////////
/* 
 * this is a ravioli class for near wall heat transfer characteristics.
 * for the moment, interfacial heat transfer resistance can be set as well
 * as a uniform dirac source, same everywhere. in the future, the source
 * should be expanded to non-uniform and the ihtr should perhaps be moved. 
 */
class HTWallModel {
  public:
    HTWallModel() : HTWallModel(0.,0.) {}

    static const HTWallModel None() {
      return HTWallModel();
    }

    static const HTWallModel Resistance(const real resist) {
      boil::oout<<"HTWallModel:near_wall_interfacial_resistance= "
                <<resist<<boil::endl;
      return HTWallModel(resist,0.);
    }

    static const HTWallModel Source(const real source) {
      boil::oout<<"HTWallModel:dirac_wall_source= "
                <<source<<boil::endl;
      return HTWallModel(0.,source);
    }

    static const HTWallModel Full(const real resist, const real source) {
      boil::oout<<"HTWallModel:near_wall_interfacial_resistance= "
                <<resist<<boil::endl;
      boil::oout<<"HTWallModel:dirac_wall_source= "
                <<source<<boil::endl;
      return HTWallModel(resist,source);
    }

    inline void set_near_wall_interfacial_resistance(const real nwir) {
      near_wall_resist = nwir;
      boil::oout<<"HTWallModel:near_wall_interfacial_resistance= "
                <<nwir<<boil::endl;
    }

    inline void set_dirac_wall_source(const real dws) {
      dirac_wall_source = dws;
      boil::oout<<"HTWallModel:dirac_wall_source= "
                <<dws<<boil::endl;
    }

    inline real temperature_node(const real len_s, const real lam_s, 
                                 const real tmp_s, const real len_f,
                                 const real lam_f, const real tmp_f) const {;
/***************************************************************************//**
*  \brief calculate temperature at node point
*             len_s         len_f
*             lam_s         lam_f
*         *-------------*------------*
*       tmp_s       tmp_node        tmp_f
*******************************************************************************/
      //return (len_f*lam_s*tmp_s + len_s*lam_f*tmp_f)/(len_f*lam_s + len_s*lam_f);
      return temperature_node(0.,len_s/lam_s, tmp_s, len_f/lam_f, tmp_f);
    }
    inline real temperature_node(const real R_s, const real tmp_s,
                                 const real R_f, const real tmp_f) const {
      //return (R_f*tmp_s+R_s*tmp_f)/(R_f+R_s);
      return temperature_node(0.,R_s,tmp_s,R_f,tmp_f);
    };
    /* with a source */
    inline real temperature_node(const real Q,
                                 const real R_s, const real tmp_s,
                                 const real R_f, const real tmp_f) const {
      return (R_f*R_s*Q+R_f*tmp_s+R_s*tmp_f)/(R_f+R_s);
    };

    real near_wall_resist;
    /* units W/m2 */
    real dirac_wall_source;

  private:

    explicit HTWallModel(const real resist, const real source) {
      near_wall_resist=resist;
      dirac_wall_source=source;
    }
};

#endif
