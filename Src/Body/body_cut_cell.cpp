#include "body.h"
#include "../Global/global_approx.h"
#include "../Field/Scalar/scalar.h"
#include "../Plot/plot.h"
//#define DEBUG

/******************************************************************************/
CutCell * Body::cut_cell(const int poly_index,
                         const int ic, const int jc, const int kc,
                         const real x0, const real y0, const real z0,
                         const real dx, const real dy, const real dz,
                         const real xc, const real yc, const real zc,
                         const real xw, const real xe,
                         const real ys, const real yn,
                         const real zb, const real zt,
                         int mflag,
                         Body * new_faces) {
/*-----------------------------------------------------------------------+
|  this routine cuts the cell and creates polygons.                      |
|         (currently used only for plotting)                             |
|                                                                        |
|  it works as follows:                                                  |
|                                                                        |
|  1. the cell is first cut with Body surface to get a cutting polygon,  |
|  2. this cutting polygon is used to create a cutting plane,            |
|  3. cutting plane is used to really cut the cell.                      |
|                                                                        | 
|  it is effectivelly a PLIC approximation.                              |
|                                                                        |
+-----------------------------------------------------------------------*/

/*----------------------------------------------------------+ 
|                                                           |
|          6-------------------7      x_01 == x_[0][0]      | 
|         /|                  /|      x_23 == x_[1][0]      |
|        / |       T         / |      x_45 == x_[0][1]      |
|       /  |       |        /  |      x_67 == x_[1][1]      |
|      4-------------------5   |                            |
|      |   |       | N     |   |      y_02 == y_[0][0]      |
|      |   |       |/      |   |      y_13 == y_[1][0]      |
|      | W---------C-------|-E |      y_46 == y_[0][1]      |
|      |   |      /|       |   |      y_57 == y_[1][1]      |
|      |   |     S |       |   |                            |
|      |   2-------|-------|---3      z_04 == z_[0][0]      |
|      |  /        |       |  /       z_15 == z_[1][0]      |
|      | /         B       | /        z_26 == z_[0][1]      |
|      |/                  |/         z_37 == z_[1][1]      |
|      0-------------------1                                |
|                                                           |
+----------------------------------------------------------*/
  const int W_E = 0;
  const int S_N = 1;
  const int B_T = 2;

#ifdef DEBUG
  std::cout<<"body_cut_cell:ic,jc,kc= "<<ic<<" "<<jc<<" "<<kc<<"\n";
#endif

/* debugging 
  bool gotcha=false;
  if( approx( z0+0.5*dz, 0.3 ) &&
      approx( y0+0.5*dy, 0.3 ) &&
      approx( x0+0.5*dz, -0.57) ) {gotcha = true; OMS(gotcha);}
*/

  /*----------------------------------------------------------------+
  |  initialize cutting cell. it is not present (NULL) by default.  |
  +----------------------------------------------------------------*/
  CutCell * ccell = NULL;

  /* new polygon coordinates (size 6 should be enough, but 
     in case cell lies on the body one might get 12 nodes) */
  real xp[12];
  real yp[12];
  real zp[12];

  /* cell coordinates */
  real x[2] = {x0, x0+dx};
  real y[2] = {y0, y0+dy};
  real z[2] = {z0, z0+dz};

  real nx=0.0;
  real ny=0.0;
  real nz=0.0;

  int n=0; /* node counter */

  Plane * p;     /* cutting plane */

  int face_match[3][2] = {0,0,0,0,0,0};  /* matching of cell face and polygon */
  real face_match_norm[3][2][3] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int nface_match=0;

  /*------------------------------------------+
  |                                           |
  |  cut cells in "i", "j" and "k" direction  |
  |                                           |
  |  generally, cells can be cut with more    |
  |   than one planar element, which may      |
  |       result in "distorted" cuts.         |
  |                                           |
  +------------------------------------------*/
  int cp;          /* cutting polygon */
  real xi, yj, zk; /* intersections */

  /* search for cuts in "i" direction */
  for(int j=0; j<2; j++)
    for(int k=0; k<2; k++)
      if(cross_seg_x( poly_index,
                      Range<real> (x[0]-tol,x[1]+tol), y[j], z[k], &xi, &cp)) {
        xp[n] = xi; yp[n] = y[j]; zp[n] = z[k];
        n++; 
        nx+=polys[cp].n(0); 
        ny+=polys[cp].n(1); 
        nz+=polys[cp].n(2);
      }

  /* search for cuts in "j" direction */
  for(int i=0; i<2; i++)
    for(int k=0; k<2; k++)
      if(cross_seg_y( poly_index,
                      x[i], Range<real> (y[0]-tol,y[1]+tol), z[k], &yj, &cp)) {
        xp[n] = x[i]; yp[n] = yj; zp[n] = z[k];
        n++; 
        nx+=polys[cp].n(0);    
        ny+=polys[cp].n(1); 
        nz+=polys[cp].n(2);
      }

  /* search for cuts in "k" direction */
  for(int i=0; i<2; i++)
    for(int j=0; j<2; j++)
      if(cross_seg_z( poly_index,
                      x[i], y[j], Range<real> (z[0]-tol,z[1]+tol), &zk, &cp)) {
        xp[n] = x[i]; yp[n] = y[j]; zp[n] = zk;
        n++; 
        nx+=polys[cp].n(0);
        ny+=polys[cp].n(1);
        nz+=polys[cp].n(2);
      }

  /* check degeneracy */
  if(n>=3){
    cut_degen( xp, yp, zp, n, ic, jc, kc);
  }

  /* cell is touched to cut surface */
  if(n<=2){
    return ccell;
  }

  /* check match-face */
  match_face( poly_index, ic, jc, kc, x, y, z, face_match, face_match_norm
           , nface_match);

  /* special treatment */
  if(nface_match>=2){
  //std::cout<<"cut_cell: nface_match=b "<<ic<<" "<<jc<<" "<<nface_match<<"\n";
    if(face_match[0][0]==1&&face_match_norm[0][0][0]<0.0) return ccell;
    if(face_match[0][1]==1&&face_match_norm[0][1][0]>0.0) return ccell;
    if(face_match[1][0]==1&&face_match_norm[1][0][1]<0.0) return ccell;
    if(face_match[1][1]==1&&face_match_norm[1][1][1]>0.0) return ccell;
    if(face_match[2][0]==1&&face_match_norm[2][0][2]<0.0) return ccell;
    if(face_match[2][1]==1&&face_match_norm[2][1][2]>0.0) return ccell;

    ccell = new CutCell();
    ccell->ijk(ic,jc,kc);

    // initialize as 1.0
    ccell->fV(1.0);
    ccell->fS(1.0,W_E,0); ccell->fS(1.0,W_E,1);
    ccell->fS(1.0,S_N,0); ccell->fS(1.0,S_N,1);
    ccell->fS(1.0,B_T,0); ccell->fS(1.0,B_T,1);
    ccell->fdxw(1.0);     ccell->fdxe(1.0);
    ccell->fdys(1.0);     ccell->fdyn(1.0);
    ccell->fdzb(1.0);     ccell->fdzt(1.0);

    for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
        for(int k=0; k<2; k++)
          ccell->fE(1.0,i,j,k);

    std::cout<<"cut_cell: nface_match=a "<<ic<<" "<<jc<<" "<<nface_match<<"\n";
    int nn=4;
    /* i-min */
    if(face_match[0][0]==1){
      xp[0]=x0;  yp[0]=y0;     zp[0]=z0;
      xp[1]=x0;  yp[1]=y0;     zp[1]=z0+dz;
      xp[2]=x0;  yp[2]=y0+dy;  zp[2]=z0+dz;
      xp[3]=x0;  yp[3]=y0+dy;  zp[3]=z0;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,W_E,0);
      ccell->fdxw(0.5);
      int i=0;
      for(int j=0; j<2; j++)
        for(int k=0; k<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    /* i-max */
    if(face_match[0][1]==1){
      xp[0]=x0+dx;  yp[0]=y0;     zp[0]=z0;
      xp[1]=x0+dx;  yp[1]=y0+dy;  zp[1]=z0;
      xp[2]=x0+dx;  yp[2]=y0+dy;  zp[2]=z0+dz;
      xp[3]=x0+dx;  yp[3]=y0;     zp[3]=z0+dz;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,W_E,1);
      ccell->fdxe(0.5);
      int i=1;
      for(int j=0; j<2; j++)
        for(int k=0; k<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    /* j-min */
    if(face_match[1][0]==1){
      xp[0]=x0;     yp[0]=y0;  zp[0]=z0;
      xp[1]=x0+dx;  yp[1]=y0;  zp[1]=z0;
      xp[2]=x0+dx;  yp[2]=y0;  zp[2]=z0+dz;
      xp[3]=x0;     yp[3]=y0;  zp[3]=z0+dz;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,S_N,0);
      ccell->fdys(0.5);
      int j=0;
      for(int i=0; i<2; j++)
        for(int k=0; k<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    /* j-max */
    if(face_match[1][1]==1){
      xp[0]=x0;     yp[0]=y0+dy;  zp[0]=z0;
      xp[1]=x0;     yp[1]=y0+dy;  zp[1]=z0+dz;
      xp[2]=x0+dx;  yp[2]=y0+dy;  zp[2]=z0+dz;
      xp[3]=x0+dx;  yp[3]=y0+dy;  zp[3]=z0;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,S_N,1);
      ccell->fdyn(0.5);
      int j=1;
      for(int i=0; i<2; j++)
        for(int k=0; k<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    /* k-min */
    if(face_match[2][0]==1){
      xp[0]=x0;     yp[0]=y0;     zp[0]=z0;
      xp[1]=x0   ;  yp[1]=y0+dy;  zp[1]=z0;
      xp[2]=x0+dx;  yp[2]=y0+dy;  zp[2]=z0;
      xp[3]=x0+dx;  yp[3]=y0;     zp[3]=z0;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,B_T,0);
      ccell->fdzb(0.5);
      int k=0;
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    if(face_match[2][1]==1){
      xp[0]=x0;     yp[0]=y0;     zp[0]=z0;
      xp[1]=x0+dx;  yp[1]=y0;     zp[1]=z0;
      xp[2]=x0+dx;  yp[2]=y0+dy;  zp[2]=z0;
      xp[3]=x0   ;  yp[3]=y0+dy;  zp[3]=z0;
      Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
      if(new_faces) new_faces->polys.push_back(*cut_1);
      ccell->fS(0.0,B_T,1);
      ccell->fdzt(0.5);
      int k=1;
      for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
          ccell->fE(0.0,i,j,k);
    }
    return ccell;
  } else if(n>=7){
    boil::aout<<"body_cut_cell:Error! (n>=7)  n= "<<n<<"\n";
    boil::aout<<"Body shape is too complicated. Use finer grid.\n";
    boil::aout<<"position of cell celter (x,y,z)= ( "<<0.5*(x[0]+x[1])<<" , "
              <<0.5*(y[0]+y[1])<<" , "<<0.5*(z[0]+z[1])<<" )\n";
    boil::aout<<"dx,dy,dz= "<<dx<<" "<<dy<<" "<<dz<<"\n";
    boil::aout<<"position (i,j,k,irank)= ( "<<ic<<" , "
              <<jc<<" , "<<kc<<" , "<<boil::cart.iam()<<" )\n";
    exit(0);
  }

  /* compute averaged normal */
  nx /= (real)n;
  ny /= (real)n;
  nz /= (real)n;

#ifdef DEBUG
  std::cout<<"body_cut_cell:1st cut n= "<<n<<"\n";
#endif

  /*-------------------------------------------------+
  |                                                  |
  |  cut cells in "i", "j" and "k" direction again,  |
  |      but with a plane created from polygon.      |
  |                                                  |
  |   that is a method to avoid "distorted" cuts.    |
  |                                                  |
  +-------------------------------------------------*/
  int nn = 0;  /* new number of nodes */

  real x_[2][2] = {-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};
  real y_[2][2] = {-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};
  real z_[2][2] = {-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX};

  /* number of times a face is cut */
  int cnt_f[3][2] = {0,0,0,0,0,0};

  /* node in fluid */
  int  n_in_fluid[2][2][2];  //-1:unknown, 0:solid, 1:fluid

  /* edge cut/non-cut */
  bool edge_cut[3][2][2];  //true:cut, false:non-cut

  /* it is quite important to initialize this */
  for(int i=0; i<2; i++){
    for(int j=0; j<2; j++){
      for(int k=0; k<2; k++){
        n_in_fluid[i][j][k] = -1;
      }
    }
  }

  for(int i=0; i<=2; i++){
    for(int j=0; j<=1; j++){
      for(int k=0; k<=1; k++){
        edge_cut[i][j][k]=false;
      }
    }
  }

  real fcx,fcy,fcz;  //face center

  if(n >= 3 && n <= 6) {

    /* cutting polygon (needed only for center of gravity) */
    arrange_poly(n, xp,yp,zp, nx,ny,nz);
    Polygon * cut_0 = new Polygon(n, xp, yp, zp); 

    /* cutting plane (uses cutting polygons center of gravity) */
    p = new Plane( cut_0->xg(), cut_0->yg(), cut_0->zg(), 
                   cut_0->n(0), cut_0->n(1), cut_0->n(2) );

    fcx=cut_0->xg();
    fcy=cut_0->yg();
    fcz=cut_0->zg();

    /* cuts in "i" direction */
    for(int j=0; j<2; j++) {
      for(int k=0; k<2; k++) {
        real xi = p->x( y[j], z[k] );
        if( (xi>x[0]-tol) && (xi<x[1]+tol) ) {
          cnt_f[S_N][j]++; 
          cnt_f[B_T][k]++; 
          x_[j][k] = xi;
          xp[nn] = xi;
          yp[nn] = y[j];
          zp[nn] = z[k];
          n_in_fluid[0][j][k] = 0;
          n_in_fluid[1][j][k] = 0;
          if ( fabs(x[0]-xi) < tol ) {
            n_in_fluid[0][j][k] = 2;
            x_[j][k] = x[0];
            xp[nn] = x[0];
          } else if ( (x[0]-xi)*p->n(0) > 0 ) {
            n_in_fluid[0][j][k] = 1;
          }
          if ( fabs(x[1]-xi) < tol ){
            n_in_fluid[1][j][k] = 2;
            x_[j][k] = x[1];
            xp[nn] = x[1];
          } else if( (x[1]-xi)*p->n(0) > 0 ) { 
            n_in_fluid[1][j][k] = 1;
          }
          if( n_in_fluid[0][j][k]+n_in_fluid[1][j][k]==1 ) {
            edge_cut[0][j][k]=true;
          }
          nn++;
        }
      }
    }

    /* cuts in "j" direction */
    for(int i=0; i<2; i++) {
      for(int k=0; k<2; k++) {
        real yj = p->y( x[i], z[k] );

        if( (yj>y[0]-tol) && (yj<y[1]+tol) ) {
          cnt_f[W_E][i]++; 
          cnt_f[B_T][k]++; 
          y_[i][k] = yj;
          xp[nn] = x[i];
          yp[nn] = yj;   
          zp[nn] = z[k];
          n_in_fluid[i][0][k] = 0;
          n_in_fluid[i][1][k] = 0;
          if ( fabs(y[0]-yj)<tol) {
            n_in_fluid[i][0][k] = 2;
            y_[i][k] = y[0];
            yp[nn] = y[0];
          } else if( (y[0]-yj)*p->n(1) > 0 ) {
            n_in_fluid[i][0][k] = 1;
          }
          if ( fabs(y[1]-yj)<tol ) {
            n_in_fluid[i][1][k] = 2;
            y_[i][k] = y[0];
            yp[nn] = y[1];
          } else if( (y[1]-yj)*p->n(1) > 0 ) {
            n_in_fluid[i][1][k] = 1;
          }
          if ( n_in_fluid[i][0][k]+n_in_fluid[i][1][k]==1 ) {
            edge_cut[1][i][k]=true;
          }
          nn++;
        }
      }
    }

    /* cuts in "k" direction */
    for(int i=0; i<2; i++) {
      for(int j=0; j<2; j++) {
        real zk = p->z( x[i], y[j] );

        if( (zk>z[0]-tol) && (zk<z[1]+tol) ) {
          cnt_f[W_E][i]++; 
          cnt_f[S_N][j]++;   
          z_[i][j] = zk;
          xp[nn] = x[i];
          yp[nn] = y[j]; 
          zp[nn] = zk;
          n_in_fluid[i][j][0] = 0;
          n_in_fluid[i][j][1] = 0;
          if ( fabs(z[0]-zk) < tol ) {
            n_in_fluid[i][j][0] = 2;
            z_[i][j] = z[0];
            zp[nn] = z[0];
          } else if( (z[0]-zk)*p->n(2) > 0 ) {
            n_in_fluid[i][j][0] = 1;
          }
          if ( fabs(z[1]-zk) < tol ) {
            n_in_fluid[i][j][1] = 2;
            z_[i][j] = z[1];
            zp[nn] = z[1];
          } else if( (z[1]-zk)*p->n(2) > 0 ) {
            n_in_fluid[i][j][1] = 1;
          }
          if ( n_in_fluid[i][j][0]+n_in_fluid[i][j][1]==1 ) {
            edge_cut[2][i][j]=true;
          }
          nn++;
        }
      }
    }
  }
  assert( nn<=12 );
  /*---------------------------------------------------+
  |  check duplication of cell-face and cutting plane  |
  +---------------------------------------------------*/
  int icflag;
  icflag = cut_duplic(p,x,y,z,n_in_fluid);
  if( icflag == 0 ){
    return ccell;
  }

  if(nn>=3){
    cut_degen( xp, yp, zp, nn, ic, jc, kc);
  }
#ifdef DEBUG
  std::cout<<"body_cut_cell:regular cut nn= "<<nn<<"\n";
#endif

  /*========================================+
  |                                         |
  |  a regular cutting face has been found  |
  |                                         |
  +========================================*/

  if(nn >= 3 && nn <= 6) {

    /*---------------------------------+
    |  cutting cell has to be created  |
    +---------------------------------*/
    ccell = new CutCell();
    ccell->ijk(ic,jc,kc);

    /*--------------------------------------------------------------------+
    |   create new cutting polygon (it is needed to compute cell volume)  |
    +--------------------------------------------------------------------*/
    //arrange_poly(n, xp,yp,zp, p->n(0),p->n(1),p->n(2));
    arrange_poly(nn, xp,yp,zp, p->n(0),p->n(1),p->n(2));
    Polygon * cut_1 = new Polygon(nn, xp, yp, zp);
#if 0
    std::cout<<"i, j, k= "<<ic<<" "<<jc<<" "<<kc<<"\n";
    for(int nnn=0; nnn<nn; nnn++){
      std::cout<<xp[nnn]<<" "<<yp[nnn]<<" "<<zp[nnn]<<"\n";
    }
#endif
    if(new_faces) new_faces->polys.push_back(*cut_1);

    /* cutting plane (uses cutting polygons center of gravity) */
    p = new Plane( cut_1->xg(), cut_1->yg(), cut_1->zg(),
                   cut_1->n(0), cut_1->n(1), cut_1->n(2) );


    /*-----------------------------------------+
    |  determine if nodal points are in fluid  |
    +-----------------------------------------*/
    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        for(int k=0; k<2; k++){
          if(n_in_fluid[i][j][k]==-1){
            real prod = p->n(0) * (x[i] - fcx)
                      + p->n(1) * (y[j] - fcy)
                      + p->n(2) * (z[k] - fcz);
            if( prod>0.0 ){
              n_in_fluid[i][j][k]=1;
            } else if ( prod<0.0 ){
              n_in_fluid[i][j][k]=0;
            } else {
              std::cout<<"body_cell_cut:Error! prod=0, irank,ic,jc,kc="<<" "
                       <<boil::cart.iam()<<" "<<ic<<" "<<jc<<" "<<kc<<"\n";
              exit(0);
            }
          }
        }
      }
    }

    for(int i=0; i<2; i++){
      for(int j=0; j<2; j++){
        for(int k=0; k<2; k++){
          if(n_in_fluid[i][j][k]==1) {
            ccell->fP(1.0,i,j,k);
          } else if(n_in_fluid[i][j][k]==0 || n_in_fluid[i][j][k]==2) {
            ccell->fP(0.0,i,j,k);
          } else {
            std::cout<<"body_cut_cell:Error! n_in_fluid)"<<"\n";
            exit(0);
          }
        }
      }
    }

    /*-------------------------+
    |  ratio of edge in fluid  |
    +-------------------------*/
    cut_edge_ratio(ccell, p, n_in_fluid, x, y, z);

    /*--------------------+
    |  face area in fluid |
    +--------------------*/
    fluid_area(ccell, n_in_fluid, edge_cut, x_, y_, z_, x, y, z, dx, dy, dz);

    /*---------------------------------------------+
    |                                              |
    |             compute cell volume              |
    |  (uses cutting polygon's center of gravity)  |
    |                                              |
    +---------------------------------------------*/
    real vol0 = dx*dy*dz;
    real vol1 = 0.0;
    /* w */ vol1 += dy * dz * ccell->fS(W_E,0) * (cut_1->xg() - x[0]);
    /* e */ vol1 += dy * dz * ccell->fS(W_E,1) * (x[1] - cut_1->xg());
    /* s */ vol1 += dx * dz * ccell->fS(S_N,0) * (cut_1->yg() - y[0]);
    /* n */ vol1 += dx * dz * ccell->fS(S_N,1) * (y[1] - cut_1->yg());
    /* b */ vol1 += dx * dy * ccell->fS(B_T,0) * (cut_1->zg() - z[0]);
    /* t */ vol1 += dx * dy * ccell->fS(B_T,1) * (z[1] - cut_1->zg());
    vol1 /= 3.0;
    real ratiov = vol1/vol0;
    real tolv = 1e-8;
    //assert( ratiov > 0.0 );
    if(ratiov<=0.0){
      std::cout<<"body_cut_cell:Error! vol1,vol0= "<<vol1<<" "<<vol0<<"\n";
      exit(0);
    }
    //assert( ratiov < 1.0 );
    if(ratiov>1.0+tolv){
      std::cout<<"body_cut_cell:Error! vol1,vol0= "<<vol1<<" "<<vol0<<" "
               <<ratiov<<"\n";
      exit(0);
    }
    if(ratiov==0.5){
      ratiov=0.5-tolv;
    }
    ccell->fV(ratiov);

    /*---------------------+
    |  check cell centers  |
    +---------------------*/
    //const real xc = 0.5 * (x[0]+x[1]); 
    //const real yc = 0.5 * (y[0]+y[1]); 
    //const real zc = 0.5 * (z[0]+z[1]); 
    
    //if( ( ( cut_1->xg() - xc ) * cut_1->n(0) +
    //      ( cut_1->yg() - yc ) * cut_1->n(1) +
    //      ( cut_1->zg() - zc ) * cut_1->n(2) ) > 0.0 ) {
      //if(ccell->fV() > 0.5) std::cout<<"body_cut_cell:fV error "
      //                               <<ccell->fV()<<"\n";
      //assert( ccell->fV() < 0.5 );
    //}

    /*-------------------------------------------+
    |  compute new (corrected) dxc, dyc and dzc  |
    +-------------------------------------------*/

    /* compute intersection. warning: it can be outside 
       the cell but i see no other way to compute it */
    const real xi = p->x(yc, zc); 
    const real yi = p->y(xc, zc); 
    const real zi = p->z(xc, yc); 

    //assert( xi != xc );
    //assert( yi != yc );
    //assert( zi != zc );
    //if( xi==xc || yi==yc || zi==zc ){
    //  std::cout<<"cut_cell:fV= "<<ccell->fV()<<" tolv= "<<tolv<<"\n";
    //  assert( fabs(ccell->fV()-0.5) < tolv);
    //}

    real fdxw=1.0, fdxe=1.0, 
         fdys=1.0, fdyn=1.0, 
         fdzb=1.0, fdzt=1.0;

    /* "i" direction; dxw & dxe */
    const real dxw = xc - xw;
    const real dxe = xe - xc;
    if( xi < xc && xi > xw ) fdxw = std::min(1.0,(xc - xi)/dxw);
    if( xi > xc && xi < xe ) fdxe = std::min(1.0,(xi - xc)/dxe);
    assert( fdxw > 0.0 );
    assert( fdxe > 0.0 );
 
    /* "j" direction; dys & dyn */
    const real dys = yc - ys;
    const real dyn = yn - yc;
    if( yi < yc && yi > ys ) fdys = std::min(1.0,(yc - yi)/dys);
    if( yi > yc && yi < yn ) fdyn = std::min(1.0,(yi - yc)/dyn);
    assert( fdys > 0.0 );
    assert( fdyn > 0.0 );

    /* "k" direction; dzb & dzt */
    const real dzb = zc - zb; 
    const real dzt = zt - zc;
    if( zi < zc && zi > zb ) fdzb = std::min(1.0,(zc - zi)/dzb);
    if( zi > zc && zi < zt ) fdzt = std::min(1.0,(zi - zc)/dzt);
    assert( fdzb > 0.0 );
    assert( fdzt > 0.0 );

    /* set the values */
    ccell->fdxw(fdxw);
    ccell->fdxe(fdxe);
    ccell->fdys(fdys);
    ccell->fdyn(fdyn);
    ccell->fdzb(fdzb);
    ccell->fdzt(fdzt);
  } /* nn>=3 && nn<=6 */

#ifdef DEBUG
  std::cout<<"body_cut_cell:end. ic,jc,kc= "<<ic<<" "<<jc<<" "<<kc<<"\n";
#endif

  return ccell;
}

/*-----------------------------------------------------------------------------+
 '$Id: body_cut_cell.cpp,v 1.34 2016/02/11 13:05:45 sato Exp $'/
+-----------------------------------------------------------------------------*/
