#include <ctime>

/*------------------------+
|  Define new MPI macros  |
+------------------------*/
#define par_comm                     int
#define par_datatype                 int
#define par_proc_null                -1
#define par_real                     0
#define par_int                      0
#define par_bool                     0
#define par_status                   int
#define par_request                  int
#define par_init(a,b)           
#define par_stop()           
#define par_wait(a,b)                
#define par_comm_rank(a,b)           *b=0
#define par_comm_size(a,b)           *b=1
#define par_barrier()
#define par_sum_real(a,b,c)          *b=*a
#define par_sum_real_n(a,b,n,c)      for(int I=0;I<(n);I++) *((b)+I)=*((a)+I)
#define par_sum_int(a,b,c)           *b=*a
#define par_sum_int_n(a,b,n,c)       for(int I=0;I<(n);I++) *((b)+I)=*((a)+I)
#define par_max_int_n(a,b,n,c)       for(int I=0;I<(n);I++) *((b)+I)=*((a)+I)
#define par_min_int_n(a,b,n,c)       for(int I=0;I<(n);I++) *((b)+I)=*((a)+I)
#define par_max_real(a,b,c)          *b=*a
#define par_min_real(a,b,c)          *b=*a
#define par_max_int(a,b,c)           *b=*a
#define par_min_int(a,b,c)           *b=*a
#define par_sendrecv(a,b,c,d,e,f,g,h,i,j,k) 
#define par_send(a,b,c,d,e) 
#define par_recv(a,b,c,d,e,f) 
#define par_isend(a,b,c,d,e,f) 
#define par_irecv(a,b,c,d,e,f)
#define par_type_vector(a,b,c,d,e)
#define par_type_commit(a)         
#define par_type_free(a)  
#define par_wtime()                  time(NULL)
