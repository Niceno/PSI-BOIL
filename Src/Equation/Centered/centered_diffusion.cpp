#include "centered.h"

/***************************************************************************//**
*  \brief Adds diffusion to right hand side.
*
*  This routine adds diffusive terms into right hand side, depending on the
*  time-stepping scheme being used. 
*
*  As an example, for Crank-Nicolson scheme, it is:
*  \f[
*      \{f\}_{old} = \{f\}_{old} + \frac{1}{2} \{D\}^{N-1}
*  \f] 
*  where \f$ \{D\} \f$ is the diffusive term and \f$ N-1 \f$ is the old time 
   step.
*
*  Diffusive term is defined as:
*  \f[
*      \{D\}^{N-1} = [A] \cdot \{ \phi \}^{N-1}
*  \f] 
*  where \f$ [A] \f$ is the system matrix containing only diffusive 
*  contributions and \f$ \{ \phi \}^{N-1} \f$ is the old vector of unknowns. 
*
*  \warning
*  Matrix \f$ [A] \f$ can change in time due to varying pysical properties.
*  To be fully consistent in such a case, this subroutine has to be invoked
*  before forming the new system of governing equations. This sequence is 
*  stipulated in main program. Obviously, these concerns are important for 
*  non-implicit time stepping schemes only.
*******************************************************************************/
void Centered::diffusion(const Scalar * diff_eddy) {

  phi.exchange();

  /* get time stepping coefficient */
  real tscn = diff_ts.N();
  real tsco = diff_ts.Nm1();
  assert( tscn > 0.0 );
  real tsc = tsco/tscn; /* 0.5 for c.n. 0.0 for fully implicit */

  /*------------------------+ 
  |  x direction (w and e)  |
  +------------------------*/
  for_ijk(i,j,k) 
    fold[i][j][k] += (tsc * A.w[i][j][k] * (phi[i-1][j][k] - phi[i][j][k]) 
                    + tsc * A.e[i][j][k] * (phi[i+1][j][k] - phi[i][j][k])); 
  
  /*------------------------+ 
  |  y direction (s and n)  |
  +------------------------*/
  for_ijk(i,j,k) 
    fold[i][j][k] += (tsc * A.s[i][j][k] * (phi[i][j-1][k] - phi[i][j][k]) 
                    + tsc * A.n[i][j][k] * (phi[i][j+1][k] - phi[i][j][k])); 

  /*------------------------+ 
  |  z direction (b and t)  |
  +------------------------*/
  for_ijk(i,j,k) 
    fold[i][j][k] += (tsc * A.b[i][j][k] * (phi[i][j][k-1] - phi[i][j][k])
                    + tsc * A.t[i][j][k] * (phi[i][j][k+1] - phi[i][j][k]));
}
