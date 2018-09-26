#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Field/Vector/vector.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
void Body::cut_stgd(const Domain & dom) {
#ifdef DEBUG
  std::cout<<"cut_stgd:start. irank= "<<boil::cart.iam()<<"\n";
#endif
 
  Body nfu,nfv,nfw; /* for checking and debugging */

  /*----------------------+
  |  cut staggered cells  |
  +----------------------*/
  for_m(m) {
#ifdef DEBUG
    std::cout<<"cut_stgd::cut m= "<<m<<"\n";
#endif
    polytags.clear();
    polytags.resize( dom.ni()*dom.nj()*dom.nk() );

    /* mark the computational cells in polygons bounding boxes */
    boil::timer.start("bounding box");
 
    /* browse through all polygons */
    for(int c=0; c<npolys(); c++) {
      /* find min and max coordinates */
      const int ifirst = vec->im(m,polys[c].minx(),tol);
      const int ilast  = vec->ip(m,polys[c].maxx(),tol);
      const int jfirst = vec->jm(m,polys[c].miny(),tol);
      const int jlast  = vec->jp(m,polys[c].maxy(),tol);
      const int kfirst = vec->km(m,polys[c].minz(),tol);
      const int klast  = vec->kp(m,polys[c].maxz(),tol);
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
    Body nf; /* for checking and debugging */
    int nccells=0;

#if 0
    std::cout<<"cut_stgd:m= "<<m<<"\n";
#endif
    for_vmijk( (*vec), m, i,j,k ) {
      const int index = i*dom.nj()*dom.nk() + j*dom.nk() + k;
      CutCell * cc = NULL;
      if( polytags[index].size() > 0 )
        cc = cut_cell( index, /* send the index of polytag structure */
                       i,j,k,
                       vec->xn(m,i),  vec->yn(m,j),  vec->zn(m,k), 
                       vec->dxc(m,i), vec->dyc(m,j), vec->dzc(m,k), 
                       vec->xc(m,i),  vec->yc(m,j),  vec->zc(m,k),
                       vec->xc(m,i-1), vec->xc(m,i+1),
                       vec->yc(m,j-1), vec->yc(m,j+1),
                       vec->zc(m,k-1), vec->zc(m,k+1),
                       0, & nf);

      if( cc ) {
        cc->ijk(i,j,k); cells[~m].push_back( *cc );
        nccells++;
      }
    }
    if(m==Comp::u()) nfu=nf;
    if(m==Comp::v()) nfv=nf;
    if(m==Comp::w()) nfw=nf;
    nccells_in[~m] = nccells;
#ifdef DEBUG
    boil::aout<<"body_cut_stgd:nccells_in= "<<nccells_in[~m]<<" "<<m<<"\n";
#endif
    boil::timer.stop("cell cutting");
  }

  /*------------+
  |  set index  | 
  +------------*/
  for_m(m) {
    for_vmijk((*vec),m,i,j,k) index[~m][i][j][k] = -1;
    for(int cc=0; cc<cells[~m].size(); cc++) {
      int i, j, k;
      cells[~m][cc].ijk(&i, &j, &k);
      index[~m][i][j][k] = cc;
    }
  }

  /*-------------+
  |  cal volume  | 
  +-------------*/
  vol_stgd(dom,nfu,nfv,nfw);

}

/*-----------------------------------------------------------------------------+
 '$Id: body_cut_stgd.cpp,v 1.2 2014/02/03 14:12:33 sato Exp $'/
+-----------------------------------------------------------------------------*/
