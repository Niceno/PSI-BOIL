/*--- related to the convective time stepping ---*/
//#define CNEW

/*--- related to diff matrix construction, no effect on accelerated discretisation ---*/
/* this is a comment from some time ago:
 * using FDM worsens the performance. Also, for FDM, the correct BndGrid for
 * Neumann boundary condition is extrapolate, rather than wall */
//#define USE_FDM_SOLID
#define USE_FDM_FLUID

/*--- related to convection old ---*/
//#define VERSION_STABLE
//#define USE_FDM_CONVECTION
