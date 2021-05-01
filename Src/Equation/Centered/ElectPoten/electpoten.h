#ifndef ELECTPOTEN_H
#define ELECTPOTEN_H

#include <cmath>
#include "../centered.h"
#include "../../../Parallel/communicator.h"
#include "../../../Parallel/Out/out.h"
//#include "../../../Global/global_outerproduct.h"
#include "../../../Global/global_vectorproduct.h"

/***************************************************************************//**
*  \brief Discretizes and solves for Poisson equation of electric potential.
*
*  The eqation is discretized in integral form:
*  \f[
*         -\int_S {\sigma_e \nabla \phi} \, dS
*        =          
*         -\int_S {\sigma_e \bf u \otimes \bf B} \, dS 
*        \; \; \; \;
*       [A]
*  \f]
*******************************************************************************/

//////////////////////////
//                      //
//  Electric Potential  //
//                      //
//////////////////////////
class ElectPoten : public Centered {
  public:
    //! Global constructor.
    /*!
        \param phi - electric potential (\f$\phi\f$),
        \param f   - right hand side array,
        \param b   - magnetic flux density,
        \param j   - electric current,
        \param v   - velocity field,
        \param v   - velocity field,
        \param t   - simulation (physical) time (\f${t}\f$),
        \param sm  - Linear solver. It acts as a solver, or as a
                     smoother for AC multigrid.
        \param mat - Holds fluid properties (\f$\sigma_e\f$).
    */
    ElectPoten(const Scalar & phi, 
             const Scalar & f,  
             Vector & b, 
             Vector & j, 
             const Vector & v, 
             Times & t, 
             Linear * sm,
             Matter * mat);
    ~ElectPoten();

    // Discretize the system of equations. 
    virtual void discretize(const Scalar * diff_eddy = NULL);

    // Computes right hand side of Poisson equation.
    real update_rhs();

    // Calculate electric current
    void update_j();

    // Output electric current
    void output_j();
    void output_j(int i);
 
    // Calculate velocity at face center
    void update_vel();

    // Calculate Lorenz force
    void force_lorenz(Vector * v);

  protected:
     
    void discretize_pressure(const Scalar * diff_eddy = NULL);
    real update_rhs_pressure();

    Vector J,B;
    Vector uf,vf,wf;

  private:
    void set_indices(Vector & v);
    void set_ranges(Vector & v, const Dir & d, const Comp & m,
                    Range<int> * i, Range<int> * j, Range<int> * k);

    void valueFace(const Vector * A,
                   const Comp m, const int i, const int j, const int k,
                   real *ax, real *ay, real *az);
};	

#endif
