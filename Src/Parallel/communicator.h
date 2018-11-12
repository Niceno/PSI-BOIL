#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#define MAX_BUFF_SIZE 1048576

#include "mpi_macros.h"
#include "tag.h"
#include "../Ravioli/dir.h"
#include "../Global/global_endl.h"
#include "../Global/global_precision.h"
#include "../Global/global_constants.h"

#include <assert.h>

/***************************************************************************//**
*  \brief Global class for controling parallel execution of the program. 
*
*  Contains subroutines for exchanging data between processors, summing up
*  the values, searching for minima, maxima and averages over processors.
*
*  \note Built upon a set of macros (see macros.h) which should make it little
*  dependent on parallel library being used. Current implementation is built 
*  upon, but analog PVM functions could also be defined inside macros.h.
*******************************************************************************/

////////////////////
//                //
//  Communicator  //
//                //
////////////////////
class Communicator {
  public:
    //! Default constructor.
    /*! Allocates memory for buffers and initializes parallel execution. */
    Communicator() {
      i = new int [MAX_BUFF_SIZE]; 
      d = new real[MAX_BUFF_SIZE];
      int     brgc=0; 
      char ** brgv=NULL;   
      par_init(&brgc, &brgv);
      par_comm_rank(MPI_COMM_WORLD, &iampro);	
      par_comm_size(MPI_COMM_WORLD, &numpro);	
      if(iam()==0) logo();
    }

    ~Communicator() {par_stop();}

    //! Acronym of "I am" returning the current processor number.
    int iam()   const {return iampro;}

    //! Returns the number of processors involved in program execution.
    int nproc() const {return numpro;}

    //! Barrier
    void barrier() const
     {par_barrier();}

    //! Exchange buffers between processors.
    void sendrecv(real * buff_send, real * buff_recv, const int size, 
                  const par_datatype & par_t, 
                  const int dir_send, const int dir_recv,
                  const Tag tag) const;
    
    //! Blocking send buffer.
    void send(real * buff_send, const int size, 
              const par_datatype & par_t, 
              const int dir_send, 
              const Tag tag) const;
    
    //! Blocking receive buffer.
    void recv(real * buff_recv, const int size, 
              const par_datatype & par_t, 
              const int dir_recv,
              const Tag tag) const;
    
    //! Non-blocking send buffer.
    void isend(real * buff_send, const int size, 
               const par_datatype & par_t, 
               const int dir_send, 
               const Tag tag,
               par_request * par) const;
    
    //! Non-blocking receive buffer.
    void irecv(real * buff_recv, const int size, 
               const par_datatype & par_t, 
               const int dir_recv,
               const Tag tag,
               par_request * req) const;
    
    //! Exchange buffers between processors.
    void sendrecv(bool * buff_send, real * buff_recv, const int size, 
                  const par_datatype & par_t, 
                  const int dir_send, const int dir_recv,
                  const Tag tag) const;
    
    //! Blocking send buffer.
    void send(bool * buff_send, const int size, 
              const par_datatype & par_t, 
              const int dir_send, 
              const Tag tag) const;
    
    //! Blocking receive buffer.
    void recv(bool * buff_recv, const int size, 
              const par_datatype & par_t, 
              const int dir_recv,
              const Tag tag) const;
    
    //! Non-blocking send buffer.
    void isend(bool * buff_send, const int size, 
               const par_datatype & par_t, 
               const int dir_send, 
               const Tag tag,
               par_request * par) const;
    
    //! Non-blocking receive buffer.
    void irecv(bool * buff_recv, const int size, 
               const par_datatype & par_t, 
               const int dir_recv,
               const Tag tag,
               par_request * req) const;

    //! Non-blocking send buffer.
    void isend(int * buff_send, const int size, 
               const par_datatype & par_t, 
               const int dir_send, 
               const Tag tag,
               par_request * par) const;
    
    //! Non-blocking receive buffer.
    void irecv(int * buff_recv, const int size, 
               const par_datatype & par_t, 
               const int dir_recv,
               const Tag tag,
               par_request * req) const;
    
    //! Send buffers from one processor to another. 
    void one2one(real * buff_send, real * buff_recv, const int size, 
                 const par_datatype & par_t, 
                 const int cpu_send, const int cpu_recv,
                 const Tag tag) const;

    //! Send buffers from one processor to another. 
    void one2one(int * buff_send, int * buff_recv, const int size,
                 const par_datatype & par_t,
                 const int cpu_send, const int cpu_recv,
                 const Tag tag) const;
    
    //! Waits for non-blocking communication to finish.
    void wait( par_request * req, par_status * stat ) const
     {par_wait( req, stat);}

    //! Summs a real number over all processors.
    int sum_real(real * a) const 
     {par_sum_real(a, d, MPI_COMM_WORLD); 
      * a = d[0]; 
      return 0;}

    //! Summs "n" real numbers over processors.
    int sum_real_n(real * a, const int n) const 
     {assert(n <= MAX_BUFF_SIZE);
      par_sum_real_n(a, d, n, MPI_COMM_WORLD); 
      for(int k=0; k<n; k++) a[k] = d[k]; 
      return 0;}

    //! Summs an integer number over all processors.
    int sum_int(int * a) const 
     {par_sum_int(a, i, MPI_COMM_WORLD); 
      * a = i[0]; 
      return 0;}

    //! Summs "n" integer numbers over processors.
    int sum_int_n(int * a, const int n) const 
     {assert(n <= MAX_BUFF_SIZE);
      par_sum_int_n(a, i, n, MPI_COMM_WORLD); 
      for(int k=0; k<n; k++) a[k] = i[k]; 
      return 0;}

    //! Returns maximum "n" integer numbers over processors.
    int max_int_n(int * a, const int n) const
     {assert(n <= MAX_BUFF_SIZE);
      par_max_int_n(a, i, n, MPI_COMM_WORLD);
      for(int k=0; k<n; k++) a[k] = i[k];
      return 0;}

    //! Returns minimum "n" integer numbers over processors.
    int min_int_n(int * a, const int n) const
     {assert(n <= MAX_BUFF_SIZE);
      par_min_int_n(a, i, n, MPI_COMM_WORLD);
      for(int k=0; k<n; k++) a[k] = i[k];
      return 0;}

    //! Returns maximum of a real number among all processors.
    int max_real(real * a) const 
     {par_max_real(a, d, MPI_COMM_WORLD); 
      * a = d[0]; 
      return 0;}

    //! Returns minimum of a real number among all processors.
    int min_real(real * a) const 
     {par_min_real(a, d, MPI_COMM_WORLD); 
      * a = d[0]; 
      return 0;}

    //! Returns maximum of a real number among all processors.
    int max_int(int * a) const 
     {par_max_int(a, i, MPI_COMM_WORLD); 
      * a = i[0]; 
      return 0;}

    //! Returns minimum of a real number among all processors.
    int min_int(int * a) const 
     {par_min_int(a, i, MPI_COMM_WORLD); 
      * a = i[0]; 
      return 0;}

    //! Returns average of a real number among all processors.
    int average_real(real * a) const 
     {par_sum_real(a, d, MPI_COMM_WORLD); 
      d[0] /= ((real) nproc());
      * a = d[0]; 
      return 0;}
    
  private:
    void logo();
    int iampro, numpro;
    int    * i;
    real * d;
};

/*----------------------+
|  global communicator  |
+----------------------*/
namespace boil {
  extern Communicator cart;
}

#endif
