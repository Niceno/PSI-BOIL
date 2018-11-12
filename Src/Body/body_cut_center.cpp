#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::cut_center(const Domain & dom) {

#ifdef DEBUG
  boil::oout<<"body_cut_center: begin\n";
#endif
 
  Body nf; /* for checking and debugging */

  /*-----------------------+
  |  cut the scalar cells  |
  +-----------------------*/
  polytags.clear();
  polytags.resize( dom.ni()*dom.nj()*dom.nk() );

  /* mark the computational cells in polygons bounding boxes */
  boil::timer.start("bounding box");

  /* browse through all polygons */
  for(int c=0; c<npolys(); c++) {
    /* find min and max coordinates */
    const int ifirst = sca->im(polys[c].minx(),tol);
    const int ilast  = sca->ip(polys[c].maxx(),tol);
    const int jfirst = sca->jm(polys[c].miny(),tol);
    const int jlast  = sca->jp(polys[c].maxy(),tol);
    const int kfirst = sca->km(polys[c].minz(),tol);
    const int klast  = sca->kp(polys[c].maxz(),tol);
    for(int i=ifirst; i<=ilast; i++)
      for(int j=jfirst; j<=jlast; j++)
        for(int k=kfirst; k<=klast; k++) {
          const int index = i*dom.nj()*dom.nk() + j*dom.nk() + k;
          polytags[index].push_back(c);
        }
  }
  boil::timer.stop("bounding box");


#ifdef DEBUG
  boil::oout<<"body_cut_center: cut cells\n";
#endif

  /* cut cells */
  boil::timer.start("cell cutting");
  int nccells=0;

  for_vijk( (*sca), i,j,k ) {
    const int index = i*dom.nj()*dom.nk() + j*dom.nk() + k;
    CutCell * cc = NULL;
    if( polytags[index].size() > 0 ){
      cc = cut_cell( index, /* send the index of polytag structure */
                     i,j,k,
                     sca->xn(i),  sca->yn(j),  sca->zn(k), 
                     sca->dxc(i), sca->dyc(j), sca->dzc(k),
                     sca->xc(i),  sca->yc(j),  sca->zc(k),
                     sca->xc(i-1), sca->xc(i+1),
                     sca->yc(j-1), sca->yc(j+1),
                     sca->zc(k-1), sca->zc(k+1),
                     1, & nf);
    }
    /* add to cut scalar cells */
    if( cc ) {
      cc->ijk(i,j,k); cells[3].push_back( *cc );
      nccells++;
    }
  }
  nccells_in[3]=nccells;
  boil::timer.stop("cell cutting");
  boil::cart.sum_int(&nccells);
  ncall_ = nccells;

  /*--------------------------------------------------------------+
  |  cut the scalar cells (buffer cells for domain decomposition) |
  +--------------------------------------------------------------*/
  polytags.clear();
  polytags.resize( dom.ni()*dom.nj()*dom.nk() );

  /* mark the computational cells in polygons bounding boxes */
  boil::timer.start("bounding box");

  /* browse through all polygons */
  for(int c=0; c<npolys(); c++) {
    /* find min and max coordinates */
    int ifirst = sca->aim(polys[c].minx(),tol);
    int ilast  = sca->aip(polys[c].maxx(),tol);
    int jfirst = sca->ajm(polys[c].miny(),tol);
    int jlast  = sca->ajp(polys[c].maxy(),tol);
    int kfirst = sca->akm(polys[c].minz(),tol);
    int klast  = sca->akp(polys[c].maxz(),tol);
    ifirst = std::max(ifirst,sifl[3]);
    ilast  = std::min(ilast ,eifl[3]);
    jfirst = std::max(jfirst,sjfl[3]);
    jlast  = std::min(jlast ,ejfl[3]);
    kfirst = std::max(kfirst,skfl[3]);
    klast  = std::min(klast ,ekfl[3]);
    for(int i=ifirst; i<=ilast; i++)
      for(int j=jfirst; j<=jlast; j++)
        for(int k=kfirst; k<=klast; k++) {
          const int index = i*dom.nj()*dom.nk() + j*dom.nk() + k;
          polytags[index].push_back(c);
        }
  }
  boil::timer.stop("bounding box");

  /* cut cells */
  boil::timer.start("cell cutting");
  nccells=0;
  for(int i=sifl[3]; i<=eifl[3]; i++){
    for(int j=sjfl[3]; j<=ejfl[3]; j++){
      for(int k=skfl[3]; k<=ekfl[3]; k++){
        if(i<(*sca).si() || i>(*sca).ei() ||
           j<(*sca).sj() || j>(*sca).ej() ||
           k<(*sca).sk() || k>(*sca).ek() ) {
          const int index = i*dom.nj()*dom.nk() + j*dom.nk() + k;

          real xw, xe, ys, yn, zb, zt;
          if (i<(*sca).si()) { xw = 2.0 * sca->xc(i) - sca->xc(i+1); }
          else               { xw = sca->xc(i-1); }
          if (i>(*sca).ei()) { xe = 2.0 * sca->xc(i) - sca->xc(i-1); }
          else               { xe = sca->xc(i+1); }
          if (j<(*sca).sj()) { ys = 2.0 * sca->yc(j) - sca->yc(j+1); }
          else               { ys = sca->yc(j-1); }
          if (j>(*sca).ej()) { yn = 2.0 * sca->yc(j) - sca->yc(j-1); }
          else               { yn = sca->yc(j+1); }
          if (k<(*sca).sk()) { zb = 2.0 * sca->zc(k) - sca->zc(k+1); }
          else               { zb = sca->zc(k-1); }
          if (k>(*sca).ek()) { zt = 2.0 * sca->zc(k) - sca->zc(k-1); }
          else               { zt = sca->zc(k+1); }

          CutCell * cc = NULL;
          if( polytags[index].size() > 0 ){
            cc = cut_cell( index, /* send the index of polytag structure */
                           i,j,k,
                           sca->xn(i),  sca->yn(j),  sca->zn(k),
                           sca->dxc(i), sca->dyc(j), sca->dzc(k),
                           sca->xc(i),  sca->yc(j),  sca->zc(k),
                           //sca->xc(i-1), sca->xc(i+1),
                           //sca->yc(j-1), sca->yc(j+1),
                           //sca->zc(k-1), sca->zc(k+1),
                           xw, xe,
                           ys, yn,
                           zb, zt,
                           1, & nf);
          }
          /* add to cut scalar cells */
          if( cc ) {
            cc->ijk(i,j,k); cells[3].push_back( *cc );
            nccells++;
          }
        }
      }
    }
  }

  nccells_bd[3]=nccells;
  boil::timer.stop("cell cutting");

  /*------------+
  |  set index  |
  +------------*/
  for_avijk((*sca),i,j,k) index[3][i][j][k] = -1;
  for(int cc=0; cc<cells[3].size(); cc++) {
    int i, j, k;
    cells[3][cc].ijk(&i, &j, &k);
    index[3][i][j][k] = cc;
  }

  // Should take into account periodic boundary condition!!!

  /*-------------+
  |  cal volume  |
  +-------------*/
  vol_center(dom,nf);

  /*--------------------------------------------------+
  |  replace original stl body with scalar cut faces  |
  +--------------------------------------------------*/
  polys = nf.polys;

#ifdef DEBUG
  boil::oout<<"body_cut_center: end\n";
#endif

}
