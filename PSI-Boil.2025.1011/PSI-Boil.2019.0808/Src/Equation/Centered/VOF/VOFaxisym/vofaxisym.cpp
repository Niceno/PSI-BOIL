#include "vofaxisym.h"

/******************************************************************************/
VOFaxisym::VOFaxisym(const Scalar & phi,
                     const Scalar & f,
                     const Scalar & kappa,
                     const Vector & u, 
                     Times & t,
                     Linear * S) :
                     //Vector * bndclr = NULL) : /* bndclr not implemented! */
  VOF(phi,f,kappa,u,t,S,NULL),
  clr( *phi.domain() ),
  axistmp ( *phi.domain() ),
  Ktmp( *phi.domain() )
{
  if(phi.domain()->is_cartesian()) {
    boil::oout<<"Warning: Initializing axisymmetric VOF on a Cartesian "
              <<"domain!"<<boil::endl;
  }
  clr = phi.shape();
  axistmp = phi.shape();
  Ktmp = phi.shape();

  set_normal_vector_method_all(NormMethod::ElviraXZ());
  //set_wall_curv_method(CurvMethod::HFnormalXZ(),Sign::neg());

  reconstruction_tolerance = 1e-4;
  reconstruction_maxiter = 5;

  /* correction */
  topo->clr = &color();

  Heaviside * htmp = new MSaxisym(&color(),NULL,&adens);
  //Heaviside * htmp = new MarchingSquares(Comp::j(),&color(),NULL,&adens);
  delete heavi;
  heavi = htmp;
}
