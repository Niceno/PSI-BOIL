#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Field/ScalarBool/scalarbool.h"
#include "../Field/VectorBool/vectorbool.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::cut_init(const Domain & dom) {
 
  /*-------------------+
  |  initialize array  |
  +-------------------*/
  sca   = new Scalar(dom); (*sca)   = 0.0; 
  bdist = new Scalar(dom); (*bdist) = 0.0;
  vecoff = new VectorBool(dom);

  real vin[3] = {0.0, 0.0, 0.0};
  vec = new Vector(dom); (*vec) = vin;

  /* set boundary condition */
  if       (dom.is_dummy(0)) {
    sca->bc().add( BndCnd( Dir::imin(), BndType::pseudo() ) );
    sca->bc().add( BndCnd( Dir::imax(), BndType::pseudo() ) );
  } else if(dom.period(0)) {
    sca->bc().add( BndCnd( Dir::imin(), BndType::periodic() ) );
    sca->bc().add( BndCnd( Dir::imax(), BndType::periodic() ) );
  } else {
    if(dom.bnd_symmetry(Dir::imin())) {
      sca->bc().add( BndCnd( Dir::imin(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::imin(), BndType::neumann() ) );
    }
    if(dom.bnd_symmetry(Dir::imax())) {
      sca->bc().add( BndCnd( Dir::imax(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::imax(), BndType::neumann() ) );
    }
  }

  if       (dom.is_dummy(1)) {
    sca->bc().add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
    sca->bc().add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
  } else if(dom.period(1)) {
    sca->bc().add( BndCnd( Dir::jmin(), BndType::periodic() ) );
    sca->bc().add( BndCnd( Dir::jmax(), BndType::periodic() ) );
  } else {
    if(dom.bnd_symmetry(Dir::jmin())) {
      sca->bc().add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::jmin(), BndType::neumann() ) );
    }
    if(dom.bnd_symmetry(Dir::jmax())) {
      sca->bc().add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::jmax(), BndType::neumann() ) );
    }
  }

  if       (dom.is_dummy(2)) {
    sca->bc().add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
    sca->bc().add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
  } else if(dom.period(2)) {
    sca->bc().add( BndCnd( Dir::kmin(), BndType::periodic() ) );
    sca->bc().add( BndCnd( Dir::kmax(), BndType::periodic() ) );
  } else {
    if(dom.bnd_symmetry(Dir::kmin())) {
      sca->bc().add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::kmin(), BndType::neumann() ) );
    }
    if(dom.bnd_symmetry(Dir::kmax())) {
      sca->bc().add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
    } else {
      sca->bc().add( BndCnd( Dir::kmax(), BndType::neumann() ) );
    }
  }

  (*bdist)=sca->shape();

  for_m(m) {
    if       (dom.is_dummy(0)) {
      vec->bc(m).add( BndCnd( Dir::imin(), BndType::pseudo() ) );
      vec->bc(m).add( BndCnd( Dir::imax(), BndType::pseudo() ) );
    } else if(dom.period(0)) {
      vec->bc(m).add( BndCnd( Dir::imin(), BndType::periodic() ) );
      vec->bc(m).add( BndCnd( Dir::imax(), BndType::periodic() ) );
    } else {
      if(dom.bnd_symmetry(Dir::imin())) {
        vec->bc(m).add( BndCnd( Dir::imin(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::imin(), BndType::neumann() ) );
      }
      if(dom.bnd_symmetry(Dir::imax())) {
        vec->bc(m).add( BndCnd( Dir::imax(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::imax(), BndType::neumann() ) );
      }
    }

    if       (dom.is_dummy(1)) {
      vec->bc(m).add( BndCnd( Dir::jmin(), BndType::pseudo() ) );
      vec->bc(m).add( BndCnd( Dir::jmax(), BndType::pseudo() ) );
    } else if(dom.period(1)) {
      vec->bc(m).add( BndCnd( Dir::jmin(), BndType::periodic() ) );
      vec->bc(m).add( BndCnd( Dir::jmax(), BndType::periodic() ) );
    } else {
      if(dom.bnd_symmetry(Dir::jmin())) {
        vec->bc(m).add( BndCnd( Dir::jmin(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::jmin(), BndType::neumann() ) );
      }
      if(dom.bnd_symmetry(Dir::jmax())) {
        vec->bc(m).add( BndCnd( Dir::jmax(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::jmax(), BndType::neumann() ) );
      }
    }

    if       (dom.is_dummy(2)) {
      vec->bc(m).add( BndCnd( Dir::kmin(), BndType::pseudo() ) );
      vec->bc(m).add( BndCnd( Dir::kmax(), BndType::pseudo() ) );
    } else if(dom.period(2)) {
      vec->bc(m).add( BndCnd( Dir::kmin(), BndType::periodic() ) );
      vec->bc(m).add( BndCnd( Dir::kmax(), BndType::periodic() ) );
    } else {
      if(dom.bnd_symmetry(Dir::kmin())) {
        vec->bc(m).add( BndCnd( Dir::kmin(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::kmin(), BndType::neumann() ) );
      }
      if(dom.bnd_symmetry(Dir::kmax())) {
        vec->bc(m).add( BndCnd( Dir::kmax(), BndType::symmetry() ) );
      } else {
        vec->bc(m).add( BndCnd( Dir::kmax(), BndType::neumann() ) );
      }
    }

  } /* for m */
#ifdef DEBUG
  std::cout<<"body_cut_init::pass B.C.set. irank= "<<boil::cart.iam()<<"\n";
#endif

  /* set range for scalar cell */
  for_m(m){
    sifl[~m]=(*vec).si(m);
    eifl[~m]=(*vec).ei(m);
    sjfl[~m]=(*vec).sj(m);
    ejfl[~m]=(*vec).ej(m);
    skfl[~m]=(*vec).sk(m);
    ekfl[~m]=(*vec).ek(m);
  }
  sifl[3]=(*sca).si();
  eifl[3]=(*sca).ei();
  sjfl[3]=(*sca).sj();
  ejfl[3]=(*sca).ej();
  skfl[3]=(*sca).sk();
  ekfl[3]=(*sca).ek();
  if( dom.dim(Comp::i()) > 1){
    if( dom.neighbour(Dir::imin()) != par_proc_null ) sifl[3]--;
    if( dom.neighbour(Dir::imax()) != par_proc_null ) eifl[3]++;
  }
  if( dom.dim(Comp::j()) > 1){
    if( dom.neighbour(Dir::jmin()) != par_proc_null ) sjfl[3]--;
    if( dom.neighbour(Dir::jmax()) != par_proc_null ) ejfl[3]++;
  }
  if( dom.dim(Comp::k()) > 1){
    if( dom.neighbour(Dir::kmin()) != par_proc_null ) skfl[3]--;
    if( dom.neighbour(Dir::kmax()) != par_proc_null ) ekfl[3]++;
  }
#ifdef DEBUG
  std::cout<<"body_cut_init:irank,sifl,siel,,,= "<<boil::cart.iam()<<" "
           <<sifl[0]<<" "<<eifl[0]<<" "<<sjfl[0]<<" "<<ejfl[0]<<"\n";
#endif

  cells.resize(4); /* u, v, w, scalar */

  index.resize(4); /* u, v, w, scalar */

  alloc3d(&index[0], vec->ni(Comp::i()),vec->nj(Comp::i()),vec->nk(Comp::i()));
  alloc3d(&index[1], vec->ni(Comp::j()),vec->nj(Comp::j()),vec->nk(Comp::j()));
  alloc3d(&index[2], vec->ni(Comp::k()),vec->nj(Comp::k()),vec->nk(Comp::k()));
  alloc3d(&index[3], sca->ni(),  sca->nj(),  sca->nk());

  /* set constant */
  dxmin = 1.0e+300;
  for_vi((*bdist),i)
    if( bdist->dxc(i)<dxmin ) dxmin = bdist->dxc(i);
  for_vj((*bdist),j)
    if( bdist->dyc(j)<dxmin ) dxmin = bdist->dyc(j);
  for_vk((*bdist),k)
    if( bdist->dzc(k)<dxmin ) dxmin = bdist->dzc(k);
  boil::cart.min_real(&dxmin);

  tol = 1.0e-8 * dxmin;
}
