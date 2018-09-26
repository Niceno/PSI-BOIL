#include "staggered.h"

/******************************************************************************/
void Staggered::set_ranges(const Dir & d, const Comp & m, 
                           Range<int> * i, Range<int> * j, Range<int> * k) {

  if( d==Dir::imin() ) {i->first(si(m)-1); i->last(si(m)-1);}
  if( d==Dir::imax() ) {i->first(ei(m)+1); i->last(ei(m)+1);}

  if( d==Dir::jmin() ) {j->first(sj(m)-1); j->last(sj(m)-1);}
  if( d==Dir::jmax() ) {j->first(ej(m)+1); j->last(ej(m)+1);}

  if( d==Dir::kmin() ) {k->first(sk(m)-1); k->last(sk(m)-1);}
  if( d==Dir::kmax() ) {k->first(ek(m)+1); k->last(ek(m)+1);}

  /*-----------------------------------------------------------------+
  |  what follows essentially corrects momentum b.c. for tangential  |
  |  velocity components (lid driven cavity). it wouln not work if   |
  |  tangential velocity components is defined only at one part of   |
  |  the boundary, which *coincides* with the domain partition.      |
  |                                                                  |
  |  or is it? i really have to check, i do not remember how it      |
  |  exaclty works or what it does :-(                               |
  +-----------------------------------------------------------------*/

  /*--------------+
  |  imin & imax  |
  +--------------*/
  if( d == Dir::imin() || d == Dir::imax() ) {
    if(m==Comp::v()) {
      if( j->first() == 1 )           j->first( u.sj(m) );
      if( j->last()  == dom->nj()-2 ) j->last ( u.ej(m) );
    }
    if(m==Comp::w()) {
      if( k->first() == 1 )           k->first( u.sk(m) );
      if( k->last()  == dom->nk()-2 ) k->last ( u.ek(m) );
    }
  }

  /*--------------+
  |  jmin & jmax  |
  +--------------*/
  if( d == Dir::jmin() || d == Dir::jmax() ) {
    if(m==Comp::u()) {
      if( i->first() == 1 )           i->first( u.si(m) );
      if( i->last()  == dom->ni()-2 ) i->last ( u.ei(m) );
    }
    if(m==Comp::w()) {
      if( k->first() == 1 )           k->first( u.sk(m) );
      if( k->last()  == dom->nk()-2 ) k->last ( u.ek(m) );
    }
  }

  /*--------------+
  |  kmin & kmax  |
  +--------------*/
  if( d == Dir::kmin() || d == Dir::kmax() ) {
    if(m==Comp::u()) {
      if( i->first() == 1 )           i->first( u.si(m) );
      if( i->last()  == dom->ni()-2 ) i->last ( u.ei(m) );
    }
    if(m==Comp::v()) {
      if( j->first() == 1 )           j->first( u.sj(m) );
      if( j->last()  == dom->nj()-2 ) j->last ( u.ej(m) );
    }
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: staggered_set_ranges.cpp,v 1.2 2009/05/09 19:06:48 niceno Exp $'/
+-----------------------------------------------------------------------------*/
