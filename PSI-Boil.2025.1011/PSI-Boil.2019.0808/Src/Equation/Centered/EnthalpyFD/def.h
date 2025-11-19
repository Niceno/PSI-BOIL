/*--- related to the convective time stepping ---*/
//#define CNEW

/*--- related to diff matrix construction ---*/
/* using FDM worsens the performance, for sure in solid.
 * Also, for FDM, the correct BndGrid for Neumann boundary condition
 * is extrapolate, rather than wall */
//#define USE_FDM_SOLID
#define USE_FDM_FLUID

/*--- related to interfacial resistance discretisation ---*/
//#define USE_FIRST_ORDER_INTRESIST

/*--- related to convection old ---*/
//#define VERSION_STABLE
//#define USE_FDM_CONVECTION
