#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define OFF    -1       /* element is switched off */
#define ON      0 

#define UNK     0       /* unknown                       */
#define GEO     1       /* geometrical (Dirichtlet) b.c. */
#define NAT     2       /* natural     (Neumann)    b.c. -> adiabatic wall */
#define OUT     3       /* natural     (Neumann)    b.c. */

#define for_nodes(n) for(int n=0; n<node.size(); n++)
#define for_elems(e) for(int e=0; e<elem.size(); e++)
#define for_sides(s) for(int s=0; s<side.size(); s++)

#endif
