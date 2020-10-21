#include "custom.h"

namespace boil {
  void setup_plane(VOF & conc, const real nnx,
                   const real nny, const real nnz,
                   const real nalp) {

    for_vijk(conc.color(),i,j,k) {
      real xpos = conc.color().xn(i);
      real ypos = conc.color().yn(j);
      real zpos = conc.color().zn(k);
      conc.nx[i][j][k] = nnx;
      conc.ny[i][j][k] = nny;
      conc.nz[i][j][k] = nnz;
      conc.nalpha[i][j][k] = nalp-nnx*xpos-nny*ypos-nnz*zpos
                           - real(nnx<0.)*nnx*conc.color().dxc(i)
                           - real(nny<0.)*nny*conc.color().dyc(j)
                           - real(nnz<0.)*nnz*conc.color().dzc(k);
      conc.nalpha[i][j][k] /= conc.color().dxc(i);
    }
    conc.forward(conc.color());

    return;
  }

  void setup_square_xz(VOF & conc, Scalar & tmp,
                       const real x0, const real z0,
                       const real sa, const real sb) {

    tmp = -3.;
    conc.color() = 0.;

    /* side 1 */
    real nnx = sa/sqrt(sa*sa+sb*sb);
    real nny = 0.;
    real nnz = sb/sqrt(sa*sa+sb*sb);
    real nalp = nnx*x0+nnz*z0;

    setup_plane(conc,nnx,nny,nnz,nalp);

    tmp += conc.color();

#if 1
    /* side 2 */
    const real x1 = x0+sb;
    const real z1 = z0-sa;
    nnx = sb/sqrt(sa*sa+sb*sb);
    nny = 0.;
    nnz = -sa/sqrt(sa*sa+sb*sb);
    nalp = nnx*x1+nnz*z1;

    setup_plane(conc,nnx,nny,nnz,nalp);

    tmp += conc.color();

    /* side 3 */
    const real x2 = x1-sa;
    const real z2 = z1-sb;
    nnx = -sa/sqrt(sa*sa+sb*sb);
    nny = 0.;
    nnz = -sb/sqrt(sa*sa+sb*sb);
    nalp = nnx*x2+nnz*z2;

    setup_plane(conc,nnx,nny,nnz,nalp);

    tmp += conc.color();

    /* side 4 */
    const real x3 = x2-sb;
    const real z3 = z2+sa;
    nnx = -sb/sqrt(sa*sa+sb*sb);
    nny = 0.;
    nnz = sa/sqrt(sa*sa+sb*sb);
    nalp = nnx*x3+nnz*z3;

    setup_plane(conc,nnx,nny,nnz,nalp);

    tmp += conc.color();
#endif

    for_vijk(tmp,i,j,k) {
      tmp[i][j][k] = std::max(tmp[i][j][k],0.0);
    }

    tmp.exchange();
    conc.color() = tmp;
    tmp = 0.;

    return;
  }
}

