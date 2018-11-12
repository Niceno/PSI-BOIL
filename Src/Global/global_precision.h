/***************************************************************************//**
*  \brief Define the precision of the program.         
*
*  Depending on the architecture and compiler (32 or 64 bit) as well as on  
*  the desired accuracy of the code, all real (non-integer) numbers as can
*  be defined as single- or double-precision.  
*
*  Therefore, a global definition for real numbers is introduced here, 
*  denoted as "real", which is an alias for "float" (single-) or "double" 
*  (double-precision).
*
*  The desired precision, also affects the MPI commands. To that end,
*  following the same approach as for real numbers, the alias "par_real"
*  has been introduced, which is an alias for "MPI_FLOAT" (single-) or
*  "MPI_DOUBLE" (double-precision).
*******************************************************************************/

/* defintions for single-precision */
// #define real     float
// #ifdef MPI
//   #define par_real MPI_FLOAT
// #endif

/* defintions for double-precision */
#define real     double
#ifdef UseMPI
  #define par_int  MPI::INT
  #define par_real MPI::DOUBLE
  #define par_bool MPI::BOOL
  // #define par_bool MPI::CHAR  //for HPC in Imperial College London
#endif
