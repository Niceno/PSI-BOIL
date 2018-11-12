// checked on: Mon Nov 12 20:52:07 CET 2007
#include "Include/psi-boil.h"

/******************************************************************************/
main(int argc, char * argv[]) {

  boil::timer.start();

  boil::plot = new PlotGMV();

  Grid1D gx(Range<real>(0,100),  8, Periodic::no());
  Grid1D gy(Range<real>(0,100),  6, Periodic::no());
  Grid1D gz(Range<real>(0,100),  4, Periodic::no());

  /*----------------------------+
  |  create distributed fields  |
  +----------------------------*/
  Domain d(gx, gy, gz);
  
  Scalar PHI(d);

  Vector uvw(d);

  /*---------------+
  |                |
  |  check scalar  |
  |                |
  +---------------*/
  boil::oout << "SCALAR" << boil::endl;
  for(int i=0; i<PHI.ni(); i++)
    for(int j=0; j<PHI.nj(); j++)
      for(int k=0; k<PHI.nk(); k++)
        PHI[i][j][k] = 8;

  for(int i=1; i<PHI.ni()-1; i++)
    for(int j=1; j<PHI.nj()-1; j++)
      for(int k=1; k<PHI.nk()-1; k++)
        PHI[i][j][k] = boil::cart.iam();

  PHI.exchange();

  /*-----------+
  |   ^ j      |
  |   |        |
  |   +--> i   |
  +-----------*/
  MPI_Barrier(MPI_COMM_WORLD);	  
  for(int m=0; m<boil::cart.nproc(); m++) {

    if(m == boil::cart.iam()) {

      /* draw separator */
      for(int k=0; k<PHI.nk(); k++) {
        boil::aout << "+---";	      
        for(int i=0; i<PHI.ni(); i++) 
          boil::aout << "--";            
        boil::aout << "---";	      
       }
      boil::aout << "+" << boil::endl;	      

      /* print values */
      for(int j=PHI.nj()-1; j>-1; j--) {
        boil::aout << "|   ";	      
        for(int k=0; k<PHI.nk(); k++) {
          for(int i=0; i<PHI.ni(); i++) 
	    boil::aout << PHI[i][j][k] << " ";      
          boil::aout << "   |   ";	      
	}
       boil::aout << m << boil::endl;	      
      }	      
    }

    MPI_Barrier(MPI_COMM_WORLD);	  
  }

  /*---------------+
  |                |
  |  check vector  |
  |                |
  +---------------*/
  boil::oout << "VECTOR" << boil::endl;
  const Comp m = Comp::w();

  for(int i=0; i<uvw.ni(m); i++)
    for(int j=0; j<uvw.nj(m); j++)
      for(int k=0; k<uvw.nk(m); k++)
        uvw[m][i][j][k] = boil::cart.iam(); //8;

  for(int i=1; i<uvw.ni(m)-1; i++)
    for(int j=1; j<uvw.nj(m)-1; j++)
      for(int k=1; k<uvw.nk(m)-1; k++)
        uvw[m][i][j][k] = boil::cart.iam();

  uvw.exchange();

  /*-----------+
  |   ^ j      |
  |   |        |
  |   +--> i   |
  +-----------*/
  MPI_Barrier(MPI_COMM_WORLD);	  
  for(int n=0; n<boil::cart.nproc(); n++) {

    if(n == boil::cart.iam()) {

      /* draw separator */
      for(int k=0; k<uvw.nk(m); k++) {
        boil::aout << "+---";	      
        for(int i=0; i<uvw.ni(m); i++) 
          boil::aout << "--";            
        boil::aout << "---";	      
       }
      boil::aout << "+" << boil::endl;	      

      /* print values */
      for(int j=uvw.nj(m)-1; j>-1; j--) {

        boil::aout << "|   ";	      
        for(int k=0; k<uvw.nk(m); k++) {
          for(int i=0; i<uvw.ni(m); i++) 
	    boil::aout << uvw[m][i][j][k] << " ";      
          boil::aout << "   |   ";	      
	  
	}
        boil::aout << n << boil::endl;	      
      }	      
    }

    MPI_Barrier(MPI_COMM_WORLD);	  
  }

  boil::timer.stop();
  boil::timer.report();

  /* used for testing only */
  boil::plot->plot(uvw, "test", 0);
}
