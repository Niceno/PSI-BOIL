#ifndef FINE_SCALAR_H
#define FINE_SCALAR_H

#include "../../../Field/Vector/vector.h"
#include "../../../Field/Scalar/scalar.h"
#include "../../../Parallel/communicator.h"

//////////////////
//              //
//  FineScalar  //
//              //
//////////////////
class FineScalar {
  public:
    FineScalar(const Scalar & PHI, const Vector & FS,
               const Scalar * NX, const Scalar * NY, const Scalar * NZ,
               const Scalar * NALPHA) :
               faceval(*FS.domain()) {

      phi = &PHI; 
      fs = &FS;
      nx = NX;  
      ny = NY;  
      nz = NZ;  
      nalpha = NALPHA;  

      /* note that bnd conditions are handled in a special manner */

      /* initialize node-based values */
      nodeval = PHI.shape();
      nodeval.allocate(PHI.ni()+1, PHI.nj()+1, PHI.nk()+1);
      nodeval.ox(1); 
      nodeval.oy(1);
      nodeval.oz(1);

      /* initialize face-based values */
      for_m(m){
        faceval(m) = FS(m).shape();
      }

      /* initialize edge-based values */
      edgex = PHI.shape();
      edgex.allocate(PHI.ni(),   PHI.nj()+1, PHI.nk()+1);
                                edgex.oy(1); edgex.oz(1);

      edgey = PHI.shape();
      edgey.allocate(PHI.ni()+1, PHI.nj(),   PHI.nk()+1);
                    edgey.ox(1);             edgey.oz(1);

      edgez = PHI.shape();
      edgez.allocate(PHI.ni()+1, PHI.nj()+1, PHI.nk()  );
                    edgez.ox(1); edgez.oy(1);
  
      phisurf = 0.5;
      }

    /* directions */
    int wsb() { return 0; }
    int wst() { return 1; }
    int wnb() { return 2; } 
    int wnt() { return 3; } 
    int esb() { return 4; } 
    int est() { return 5; } 
    int enb() { return 6; } 
    int ent() { return 7; } 

    int ws() { return 8; }
    int wn() { return 9; }
    int wb() { return 10; }
    int wt() { return 11; }
    int es() { return 12; }
    int en() { return 13; }
    int eb() { return 14; }
    int et() { return 15; }
    int sb() { return 16; }
    int st() { return 17; }
    int nb() { return 18; }
    int nt() { return 19; }

    int w() { return 20; }
    int e() { return 21; }
    int s() { return 22; }
    int n() { return 23; }
    int b() { return 24; }
    int t() { return 25; }

    real & value(int i, int j, int k, int dir);
    void evaluate();

   private:
    typedef struct {
      bool rl; /* does the interface cross bw point and cell centre? */
      int marker; /* -1 below intface (phi=0), 0 at intface (phi=0.5) */
      real dist;
    } point_coord;

    void eval_node();
    void eval_face();
    void eval_edge();

    real edge_eval_val(const real phi1,
                       const int i1, const int j1, const int k1,
                       const int dir1,
                       const real phi2,
                       const int i2, const int j2, const int k2,
                       const int dir2);

    point_coord edge_eval_point(const int i, const int j, 
                                const int k, const int dir);

    void bdcond_node();
    void bdcond_face();
    void bdcond_x();
    void bdcond_y();
    void bdcond_z();

    Scalar nodeval; /* nodal values */
    Vector faceval; /* face values  */
    Scalar edgex; /* edge values, staggered in y and z */
    Scalar edgey; /* edge values, staggered in x and z */
    Scalar edgez; /* edge values, staggered in x and y */

    const Scalar * phi;
    const Vector * fs;
    const Scalar * nx;
    const Scalar * ny;
    const Scalar * nz;
    const Scalar * nalpha;
 
    real phisurf;
};

#endif
