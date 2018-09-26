#include "mpi.h"
#include <ctime>

#define UseMPI

/*------------------------+
|  Define new MPI macros  |
+------------------------*/
#define par_comm                     MPI_Comm
#define par_datatype                 MPI_Datatype
#define par_proc_null                MPI_PROC_NULL
#define par_status                   MPI_Status
#define par_request                  MPI_Request
#define par_init(a,b)                MPI_Init       ((a),(b))   
#define par_stop()                   MPI_Finalize   ()   
#define par_wait(a,b)                MPI_Wait       ((a),(b))   
#define par_comm_rank(a,b)           MPI_Comm_rank  ((a),(b))   
#define par_comm_size(a,b)           MPI_Comm_size  ((a),(b))   
#define par_barrier()                MPI_Barrier    (MPI_COMM_WORLD)
#define par_sum_real(a,b,c)     \
        MPI_Allreduce((a),(b),1,par_real,MPI_SUM,(c))
#define par_sum_real_n(a,b,n,c) \
        MPI_Allreduce((a),(b),(n),par_real,MPI_SUM,(c))
#define par_sum_int(a,b,c)        \
        MPI_Allreduce((a),(b),1,par_int,MPI_SUM,(c))
#define par_sum_int_n(a,b,n,c)    \
        MPI_Allreduce((a),(b),(n),par_int,MPI_SUM,(c))
#define par_max_int_n(a,b,n,c)    \
        MPI_Allreduce((a),(b),(n),par_int,MPI_MAX,(c))
#define par_min_int_n(a,b,n,c)    \
        MPI_Allreduce((a),(b),(n),par_int,MPI_MIN,(c))
#define par_max_real(a,b,c)     \
        MPI_Allreduce((a),(b),1,par_real,MPI_MAX,(c))
#define par_min_real(a,b,c)     \
        MPI_Allreduce((a),(b),1,par_real,MPI_MIN,(c))
#define par_max_int(a,b,c)     \
        MPI_Allreduce((a),(b),1,par_int,MPI_MAX,(c))
#define par_min_int(a,b,c)     \
        MPI_Allreduce((a),(b),1,par_int,MPI_MIN,(c))
#define par_sendrecv(a,b,c,d,e,f,g,h,i,j,k) \
        MPI_Sendrecv((a),(b),(c),(d),(e),(f),(g),(h),(i),(j), \
                     MPI_COMM_WORLD,(k)) 
#define par_send(a,b,c,d,e) \
        MPI_Send((a),(b),(c),(d),(e),MPI_COMM_WORLD) 
#define par_recv(a,b,c,d,e,f) \
        MPI_Recv((a),(b),(c),(d),(e),MPI_COMM_WORLD,(f)) 
#define par_isend(a,b,c,d,e,f) \
        MPI_Isend((a),(b),(c),(d),(e),MPI_COMM_WORLD,(f)) 
#define par_irecv(a,b,c,d,e,f) \
        MPI_Irecv((a),(b),(c),(d),(e),MPI_COMM_WORLD,(f)) 
#define par_type_vector(a,b,c,d,e)   MPI_Type_vector((a),(b),(c),(d),(e))
#define par_type_commit(a)           MPI_Type_commit((a))
#define par_type_free(a)             MPI_Type_free  ((a))
#define par_wtime()                  MPI_Wtime      ()

/*-----------------------------------------------------------------------------+
 '$Id: parallel.h,v 1.10 2015/05/05 14:11:44 sato Exp $'/
+-----------------------------------------------------------------------------*/
