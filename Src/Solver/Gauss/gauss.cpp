#include "gauss.h"

#define POS(i,j,k) ( (i-1)*nj*nk + (j-1)*nk + (k-1) )

/******************************************************************************/
void Gauss::solve(Matrix & A, Scalar & X, Scalar & B) {

  real ** a;
  real  * x;
  real  * b;
  int   * indx;

  const int ni = X.ni()-2;
  const int nj = X.nj()-2;
  const int nk = X.nk()-2;

  int n = ni*nj*nk;

  /*------------------------------------------+ 
  |  get number of cells on other processors  |
  +------------------------------------------*/
  int * nd, * nf;
  alloc1d(&nd, boil::cart.nproc()); // number of cells in each domain
  alloc1d(&nf, boil::cart.nproc()); // first cell in each domain (counting) */
  nd[ boil::cart.iam() ] = n;

  boil::cart.sum_int_n(nd, boil::cart.nproc());

  /*----------------------------------------------+
  |  compute number of cells over all processors  |
  |      and first cell for each processor        |
  +----------------------------------------------*/
  n       = 0; for(int i=0; i<boil::cart.nproc(); i++) n += nd[i];
  nf[ 0 ] = 0; for(int i=1; i<boil::cart.nproc(); i++) nf[i] = nf[i-1]+nd[i-1];

  /*-------------------------------+
  |  look for indices at periodic  |
  |    and processor boundaries    |
  +-------------------------------*/
  Scalar pos(*X.domain()); 

  /* set the proper numbers inside the domain and -1's on the boundary */
  for_avijk(X,i,j,k) pos[i][j][k] = -1.01;
  for_vijk(X,i,j,k)  pos[i][j][k] = POS(i,j,k) + 0.01 + nf[ boil::cart.iam() ];

  /* get indices from other processors */
  pos.exchange();

  /*--------------------------+
  |  allocate memory for all  |
  +--------------------------*/
  alloc2d(&a, n, n);
  alloc1d(&b, n);
  alloc1d(&x, n);
  alloc1d(&indx ,n);

  int c, cp;

  /*-------------------+
  |  form full matrix  |
  +-------------------*/
  for_vijk(X, i,j,k) {

    c = (int)pos[i][j][k];

    a[c][c] = A.c[i][j][k];
    
    /* w */
    cp = (int)pos[i-1][j][k];
    if(cp > -1) a[c][cp] = -A.w[i][j][k];
    else        b[c]     += A.w[i][j][k] * X[i-1][j][k];

    /* e */
    cp = (int)pos[i+1][j][k];
    if(cp > -1) a[c][cp] = -A.e[i][j][k];
    else        b[c]     += A.e[i][j][k] * X[i+1][j][k];

    /* s */
    cp = (int)pos[i][j-1][k];
    if(cp > -1) a[c][cp] = -A.s[i][j][k];
    else        b[c]     += A.s[i][j][k] * X[i][j-1][k];

    /* n */
    cp = (int)pos[i][j+1][k];
    if(cp > -1) a[c][cp] = -A.n[i][j][k];
    else        b[c]     += A.n[i][j][k] * X[i][j+1][k];

    /* b */
    cp = (int)pos[i][j][k-1];
    if(cp > -1) a[c][cp] = -A.b[i][j][k];
    else        b[c]     += A.b[i][j][k] * X[i][j][k-1];

    /* t */
    cp = (int)pos[i][j][k+1];
    if(cp > -1) a[c][cp] = -A.t[i][j][k];
    else        b[c]     += A.t[i][j][k] * X[i][j][k+1];
  }

  /*-------------------+
  |  form vectors too  |
  +-------------------*/
  for_vijk(X, i,j,k) 
    b[(int)pos[i][j][k]] += B[i][j][k];

  // plot_system(a,b,n,nj,nk,"A-separate.eps");

  /*-----------------------+
  |  summ the matrices up  |
  +-----------------------*/
  boil::cart.sum_real_n(&b[0],n);
  boil::cart.sum_real_n(&a[0][0],n*n);

  // plot_system(a,b,n,nj,nk,"A-summedup.eps");

  /*-----------+
  |  solve it  |
  +-----------*/
  legs (a,b,x,indx,n);

  // plot_system(a,b,n,nj,nk,"A-after.eps");

  /*-------------------------------------+
  |  copy solution back to input vector  |
  +-------------------------------------*/
  for_vijk(X, i,j,k)
    X[i][j][k] = x[(int)pos[i][j][k]];

  /*--------------+
  |  free memory  |
  +--------------*/
  dealloc2d(&a);
  dealloc1d(&x);
  dealloc1d(&b);
  dealloc1d(&indx);

  //pos.deallocate();
  delete [] nd;
  delete [] nf;
}

/****************************************************************************/
void Gauss::legs (real ** a, real * b, real * x, int * indx, int n) {
/*----------------------------------------------------------+
|  Function to solve the equation a[][] x[] = b[] with the  |
|  partial-pivoting Gaussian elimination scheme.            |
|  Copyright (c) Tao Pang 2001.                             |
+----------------------------------------------------------*/

  int i,j;

  elgs (a,indx,n);

  for(i = 0; i < n-1; ++i)
  {
    for(j = i+1; j < n; ++j)
    {
      b[indx[j]] = b[indx[j]]-a[indx[j]][i]*b[indx[i]];
    }
  }

  x[n-1] = b[indx[n-1]]/a[indx[n-1]][n-1];
  for (i = n-2; i>=0; i--)
  {
    x[i] = b[indx[i]];
    for (j = i+1; j < n; ++j)
    {
      x[i] = x[i]-a[indx[i]][j]*x[j];
    }
    x[i] = x[i]/a[indx[i]][i];
  }
}

/****************************************************************************/
void Gauss::elgs(real ** a, int * indx, int n) {
/*-----------------------------------------------------------------+
|  Function to perform the partial pivoting Gaussian elimination.  |
|  a[][] is the original matrix in the input and transformed       |
|  matrix plus the pivoting element ratios below the diagonal      |
|  in the output.  indx[] records the pivoting order.              |
|  Copyright (c) Tao Pang 2001.                                    |
+-----------------------------------------------------------------*/

  int i, j, k, itmp;
  real c1, pi, pi1, pj;
  real * c;

  alloc1d(&c, n);

/* Initialize the index */

  for (i = 0; i < n; ++i)
  {
    indx[i] = i;
  }

/* Find the rescaling factors, one from each row */
 
  for (i = 0; i < n; ++i)
  {
    c1 = 0;
    for (j = 0; j < n; ++j)
    {
      if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
    }
    c[i] = c1;
  }

/* Search the pivoting (largest) element from each column */ 

  for (j = 0; j < n-1; ++j)
  {
    pi1 = 0;
    for (i = j; i < n; ++i)
    {
      pi = fabs(a[indx[i]][j])/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }

/* Interchange the rows via indx[] to record pivoting order */

    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < n; ++i)
    {
      pj = a[indx[i]][j]/a[indx[j]][j];

/* Record pivoting ratios below the diagonal */

      a[indx[i]][j] = pj;

/* Modify other elements accordingly */

      for (k = j+1; k < n; ++k)
      {
        a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
      }
    }
  }
  dealloc1d(&c);
}

/****************************************************************************/
void Gauss::plot_system(real ** a, real * b, int c, int nj, int nk, 
                        char * nam) {

  Board board;
  const int box= 4;
  const int dis= 8;

  /*------------------------------------------------------------+ 
  |  draw divisions by "nk" (n-s neighbours) - checkboard like  |
  +------------------------------------------------------------*/
  board.setLineWidth( 0.0 );
  board.setFillColorRGBi(224,224,224);
  board.setPenColorRGBi(224,224,224);
  for(int i=0; i<c; i += nk) {
    for(int j=0; j<c; j += nk) 
      if( i%(2*nk) - j%(2*nk) == 0 )
        board.drawRectangle( j*box, (c-i)*box, nk*box, nk*box ); /* matrix */

    if( i%(2*nk) == 0 )
      board.drawRectangle( c*box + dis, (c-i)*box, box, nk*box ); /* rhs */
  }

  /*----------------------------------------------------------------+ 
  |  draw divisions by "nk*nj" (e-w neighbours) - coarse grid like  |
  +----------------------------------------------------------------*/
  board.setLineWidth( 0.25 );
  board.setPenColorRGBi( 0, 0, 0 );
  for(int i=0; i<=c; i += nk*nj) {
    board.drawLine(0, (c-i)*box, c*box, (c-i)*box); /* matrix */
    board.drawLine(i*box, c*box, i*box, 0);         /* matrix */

    board.drawLine(c*box + dis, (c-i)*box, (c+1)*box+dis, (c-i)*box); /* rhs */
  }
  board.drawLine( c   *box + dis, c*box,  c   *box + dis, 0); /* rhs */
  board.drawLine((c+1)*box + dis, c*box, (c+1)*box + dis, 0); /* rhs */

  /*----------------------+ 
  |  draw matrix entries  |
  +----------------------*/
  board.setLineWidth( 0.25 );
  board.setPenColorRGBi( 0, 0, 0 );
  for(int i=0; i<c; i++) { // rows
    for(int j=0; j<c; j++) {
      if( a[i][j] > 0.0 ) { // colums 
        board.setFillColorRGBi(255,  0,  0);
        board.drawRectangle( j*box, c*box - i*box, box, box );
      }
      if( a[i][j] < 0.0 ) {
        board.setFillColorRGBi(  0,  0,255);
        board.drawRectangle( j*box, c*box - i*box, box, box );
      }
    }
  }

  /*-------------------+ 
  |  draw rhs entries  |
  +-------------------*/
  for(int i=0; i<c; i++) { // rows
    if( b[i] > 0.0 ) {
      board.setFillColorRGBi(255,  0,  0);
      board.drawRectangle( c*box + dis, c*box - i*box, box, box );
    }
    if( b[i] < 0.0 ) {
      board.setFillColorRGBi(  0,  0,255);
      board.drawRectangle( c*box + dis, c*box - i*box, box, box );
    }
  }

  /* file name extension */
  std::string        name(nam);
  std::ostringstream numb;
  numb << "_";
  numb.fill('0');
  numb.width(3);
  numb << boil::cart.iam();
  numb << ".eps";
  name += numb.str();

  board.saveEPS( name.c_str() );
}

/*-----------------------------------------------------------------------------+
 '$Id: gauss.cpp,v 1.9 2008/11/17 19:23:24 niceno Exp $'/
+-----------------------------------------------------------------------------*/
