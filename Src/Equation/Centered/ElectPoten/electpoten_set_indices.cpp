#include "electpoten.h"

/******************************************************************************/
void ElectPoten::set_indices(Vector & u) {

  /*-------------------------------------+ 
  |  set start and end indices properly  |
  +-------------------------------------*/
  Comp m = Comp::u();
  if( !u.bc(m).type( Dir::imin(), BndType::periodic() ) && dom->coord(m)==0)             {u.si(m)++;}
  if( !u.bc(m).type( Dir::imax(), BndType::periodic() ) && dom->coord(m)==dom->dim(m)-1) {u.ei(m)--;}

  m = Comp::v();
  if( !u.bc(m).type( Dir::jmin(), BndType::periodic() ) && dom->coord(m)==0)             {u.sj(m)++;}
  if( !u.bc(m).type( Dir::jmax(), BndType::periodic() ) && dom->coord(m)==dom->dim(m)-1) {u.ej(m)--;}

  m = Comp::w();
  if( !u.bc(m).type( Dir::kmin(), BndType::periodic() ) && dom->coord(m)==0)             {u.sk(m)++;}
  if( !u.bc(m).type( Dir::kmax(), BndType::periodic() ) && dom->coord(m)==dom->dim(m)-1) {u.ek(m)--;}

  /*-----------------------------------------------------+
  |  transform b.c. ranges from scalar to vector values  |
  +-----------------------------------------------------*/
  for_m(m)
    for( int b=0; b<u.bc(m).count(); b++ ) {

      // needed??? if( u.bc(m).exists(b) ) {

      /* these ranges are the same as for scalar ... */
      Range<int> ir(  u.bc(m).at(b).si(),  u.bc(m).at(b).ei()  ); 
      Range<int> jr(  u.bc(m).at(b).sj(),  u.bc(m).at(b).ej()  );
      Range<int> kr(  u.bc(m).at(b).sk(),  u.bc(m).at(b).ek()  );

      Dir d = u.bc(m).direction(b);
      set_ranges(u, d, m, &ir, &jr, &kr);

      u.bc(m).at(b).si( ir.first() );
      u.bc(m).at(b).sj( jr.first() );
      u.bc(m).at(b).sk( kr.first() );

      u.bc(m).at(b).ei( ir.last() );
      u.bc(m).at(b).ej( jr.last() );
      u.bc(m).at(b).ek( kr.last() );

      // needed??? }
    }

}
