#include "bndcnd.h"

/******************************************************************************/
BndCnd::BndCnd( const Domain & d ) 
 : typ( BndType::dirichlet() ),
   dir( Dir::undefined() ),
   dom( &d ),
   ir(0,-1), jr(0,-1), kr(0,-1) {

  val[0]=0.0; 
  val[1]=0.0; 
  val[2]=0.0;
}

/******************************************************************************/
BndCnd::BndCnd(const Dir & dr, const BndType & bt, const BndVal & u, 
                                                   const BndVal & v, 
                                                   const BndVal & w) 
/*----------------------------------------------------------------+
|  creates a temprorary object when sending a parameter to "add"  |
+----------------------------------------------------------------*/
 : dir(dr), 
   typ(bt), 
   ir(0,-1), jr(0,-1), kr(0,-1) {
  val[0] = u; 
  val[1] = v; 
  val[2] = w;
}

/******************************************************************************/
BndCnd::BndCnd( const Dir        & dr, 
                const Range<int> & j_r, 
                const Range<int> & k_r, const BndType & bt, const BndVal & u,
                                                            const BndVal & v, 
                                                            const BndVal & w ) 
/*----------------------------------------------------------------+
|  creates a temprorary object when sending a parameter to "add"  |
+----------------------------------------------------------------*/
 : dir(dr), 
   typ(bt),
   ir(0,-1), jr(j_r), kr(k_r) {
 
  assert( dr == Dir::imin() || dr == Dir::imax() );

  val[0] = u; 
  val[1] = v; 
  val[2] = w;
}

/******************************************************************************/
BndCnd::BndCnd( const Range<int> & i_r,
                const Dir        & dr,  
                const Range<int> & k_r, const BndType & bt, const BndVal & u,
                                                            const BndVal & v, 
                                                            const BndVal & w ) 
/*----------------------------------------------------------------+
|  creates a temprorary object when sending a parameter to "add"  |
+----------------------------------------------------------------*/
 : dir(dr), 
   typ(bt),
   ir(i_r), jr(0,-1), kr(k_r) {
 
  assert( dr == Dir::jmin() || dr == Dir::jmax() );

  val[0] = u; 
  val[1] = v; 
  val[2] = w;
}

/******************************************************************************/
BndCnd::BndCnd( const Range<int> & i_r,
                const Range<int> & j_r, 
                const Dir        & dr,  const BndType & bt, const BndVal & u,
                                                            const BndVal & v, 
                                                            const BndVal & w ) 
/*----------------------------------------------------------------+
|  creates a temprorary object when sending a parameter to "add"  |
+----------------------------------------------------------------*/
 : dir(dr), 
   typ(bt),
   ir(i_r), jr(j_r), kr(0,-1) {
 
  assert( dr == Dir::kmin() || dr == Dir::kmax() );

  val[0] = u; 
  val[1] = v; 
  val[2] = w;
}

/******************************************************************************/
BndCnd::BndCnd( const Range<int> & i_r,
                const Range<int> & j_r, 
                const Range<int> & k_r, 
                const Dir        & dr,  const BndType & bt, const BndVal & u,
                                                            const BndVal & v, 
                                                            const BndVal & w ) 
/*----------------------------------------------------------------+
|  creates a temprorary object when sending a parameter to "add"  |
+----------------------------------------------------------------*/
 : dir(dr), 
   typ(bt),
   ir(i_r), jr(j_r), kr(k_r) {
 
  assert( dr == Dir::ibody() );

  val[0] = u; 
  val[1] = v; 
  val[2] = w;
}

/******************************************************************************/
void BndCnd::modify(BndCnd bc) {
  for(int b=0; b<count(); b++){
    if(section.at(b).typ == bc.typ && section.at(b).dir == bc.dir){
      section.at(b).val[0] = bc.val[0];
      section.at(b).val[1] = bc.val[1];
      section.at(b).val[2] = bc.val[2];
    }
  }
}

/******************************************************************************/
void BndCnd::add(BndCnd bc) {

  bool fail;

  /*----------------------------------------------+
  |  1. check consistency of boundary conditions  | 
  +----------------------------------------------*/
  fail = false;
  if( bc.typ == BndType::periodic() ) {
    if( bc.dir == Dir::imin() && !dom->period(0) ) fail = true;
    if( bc.dir == Dir::imax() && !dom->period(0) ) fail = true;
    if( bc.dir == Dir::jmin() && !dom->period(1) ) fail = true;
    if( bc.dir == Dir::jmax() && !dom->period(1) ) fail = true;
    if( bc.dir == Dir::kmin() && !dom->period(2) ) fail = true;
    if( bc.dir == Dir::kmax() && !dom->period(2) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define periodic b.c. " 
               << "to a non-periodic domain. Exiting!" << boil::endl;
    exit(0);
  }

  fail = false;
  if( bc.typ != BndType::periodic() && bc.typ != BndType::pseudo() ) {
    if( bc.dir == Dir::imin() && dom->period(0) ) fail = true;
    if( bc.dir == Dir::imax() && dom->period(0) ) fail = true;
    if( bc.dir == Dir::jmin() && dom->period(1) ) fail = true;
    if( bc.dir == Dir::jmax() && dom->period(1) ) fail = true;
    if( bc.dir == Dir::kmin() && dom->period(2) ) fail = true;
    if( bc.dir == Dir::kmax() && dom->period(2) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define real non-periodic b.c. " 
               << "to a periodic domain. Exiting!" << boil::endl;
    exit(0);
  }

#if 1
  fail = false;
  if( bc.typ == BndType::pseudo() ) {
    if( bc.dir == Dir::imin() && !dom->is_dummy(0) ) fail = true;
    if( bc.dir == Dir::imax() && !dom->is_dummy(0) ) fail = true;
    if( bc.dir == Dir::jmin() && !dom->is_dummy(1) ) fail = true;
    if( bc.dir == Dir::jmax() && !dom->is_dummy(1) ) fail = true;
    if( bc.dir == Dir::kmin() && !dom->is_dummy(2) ) fail = true;
    if( bc.dir == Dir::kmax() && !dom->is_dummy(2) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define a dummy b.c. "
               << "in a real direction. Exiting!" << boil::endl;
    exit(0);
  }

  fail = false;
  if( bc.typ != BndType::pseudo() ) {
    if( bc.dir == Dir::imin() && dom->is_dummy(0) ) fail = true;
    if( bc.dir == Dir::imax() && dom->is_dummy(0) ) fail = true;
    if( bc.dir == Dir::jmin() && dom->is_dummy(1) ) fail = true;
    if( bc.dir == Dir::jmax() && dom->is_dummy(1) ) fail = true;
    if( bc.dir == Dir::kmin() && dom->is_dummy(2) ) fail = true;
    if( bc.dir == Dir::kmax() && dom->is_dummy(2) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define a real b.c. "
               << "in a pseudo-direction. Exiting!" << boil::endl;
    exit(0);
  }
#endif

  fail = false;
  if( bc.typ == BndType::symmetry() ) {
    if( bc.dir == Dir::imin() && !dom->bnd_symmetry(Dir::imin()) ) fail = true;
    if( bc.dir == Dir::imax() && !dom->bnd_symmetry(Dir::imax()) ) fail = true;
    if( bc.dir == Dir::jmin() && !dom->bnd_symmetry(Dir::jmin()) ) fail = true;
    if( bc.dir == Dir::jmax() && !dom->bnd_symmetry(Dir::jmax()) ) fail = true;
    if( bc.dir == Dir::kmin() && !dom->bnd_symmetry(Dir::kmin()) ) fail = true;
    if( bc.dir == Dir::kmax() && !dom->bnd_symmetry(Dir::kmax()) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define symmetry b.c. "
               << "to a non-symmetric domain. Exiting!" << boil::endl;
    exit(0);
  }

  fail = false;
  if( bc.typ != BndType::symmetry() ) {
    if( bc.dir == Dir::imin() && dom->bnd_symmetry(Dir::imin()) ) fail = true;
    if( bc.dir == Dir::imax() && dom->bnd_symmetry(Dir::imax()) ) fail = true;
    if( bc.dir == Dir::jmin() && dom->bnd_symmetry(Dir::jmin()) ) fail = true;
    if( bc.dir == Dir::jmax() && dom->bnd_symmetry(Dir::jmax()) ) fail = true;
    if( bc.dir == Dir::kmin() && dom->bnd_symmetry(Dir::kmin()) ) fail = true;
    if( bc.dir == Dir::kmax() && dom->bnd_symmetry(Dir::kmax()) ) fail = true;
  }
  if( fail ) {
    boil::oout << "Fatal: trying to define non-symmetric b.c. "
               << "to a symmetric domain. Exiting!" << boil::endl;
    exit(0);
  }

  /*------------------------------------------------------------+
  |  1. skip this b.c. if it does not belong to this subdomain  | 
  +------------------------------------------------------------*/
  fail = false; 

  if(bc.dir == Dir::imin() && dom->coord(Comp::i()) != 0)             
    fail = true;
  if(bc.dir == Dir::imax() && dom->coord(Comp::i()) != dom->dim(Comp::i())-1) 
    fail = true;

  if(bc.dir == Dir::jmin() && dom->coord(Comp::j()) != 0)             
    fail = true;
  if(bc.dir == Dir::jmax() && dom->coord(Comp::j()) != dom->dim(Comp::j())-1) 
    fail = true;

  if(bc.dir == Dir::kmin() && dom->coord(Comp::k()) != 0)             
    fail = true;
  if(bc.dir == Dir::kmax() && dom->coord(Comp::k()) != dom->dim(Comp::k())-1) 
    fail = true;

  if( fail ) {
    bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
    bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    section.push_back(bc);
    return;
  }

#if 0
  /*---------------------------------------------------+
  |  2. set range for directions explicitelly defined  | 
  +---------------------------------------------------*/
  if(bc.dir == Dir::imin()) {bc.ir.first(boil::BW-1);
                             bc.ir.last (boil::BW-1);}
  if(bc.dir == Dir::imax()) {bc.ir.first(dom->ni()-boil::BW);
                             bc.ir.last (dom->ni()-boil::BW);}

  if(bc.dir == Dir::jmin()) {bc.jr.first(boil::BW-1);
                             bc.jr.last (boil::BW-1);}
  if(bc.dir == Dir::jmax()) {bc.jr.first(dom->nj()-boil::BW);
                             bc.jr.last (dom->nj()-boil::BW);}

  if(bc.dir == Dir::kmin()) {bc.kr.first(boil::BW-1);
                             bc.kr.last (boil::BW-1);}
  if(bc.dir == Dir::kmax()) {bc.kr.first(dom->nk()-boil::BW);
                             bc.kr.last (dom->nk()-boil::BW);}

  boil::oout<<"Bnd after 2 "<<bc.dir<<" "<<bc.ir.first()<<" "<<bc.ir.last()<<" | "<<bc.kr.first()<<" "<<bc.kr.last()<<" || "<<boil::endl;

  /*--------------------------------------+
  |  3. set the range for the directions  |
  |         not explicitelly defined      |
  +--------------------------------------*/
  Range<int> cxg = dom->cxg(); 
  Range<int> cyg = dom->cyg(); 
  Range<int> czg = dom->czg(); 

  // What is this ??? (Yohei)
  if( bc.ir.first() == 0 && bc.ir.last() == -1 ) {
    bc.ir.first( cxg.first() + boil::BW - 1 ); // as sx()
    bc.ir.last ( cxg.last()  + boil::BW - 1 ); // as ex()
  }
  if( bc.jr.first() == 0 && bc.jr.last() == -1 ) {
    bc.jr.first( cyg.first() + boil::BW - 1 ); // as sy()
    bc.jr.last ( cyg.last()  + boil::BW - 1 ); // as ey()
  }
  if( bc.kr.first() == 0 && bc.kr.last() == -1 ) {
    bc.kr.first( czg.first() + boil::BW - 1 ); // as sz()
    bc.kr.last ( czg.last()  + boil::BW - 1 ); // as ez()

  }

  boil::oout<<"Bnd after 3 "<<bc.dir<<" "<<bc.ir.first()<<" "<<bc.ir.last()<<" | "<<bc.kr.first()<<" "<<bc.kr.last()<<" || "<<boil::endl;

  /*---------------------------------------------------+
  |  4. correct the range for directions defined by 3  |
  +---------------------------------------------------*/
  bool i = bc.dir==Dir::imin() || bc.dir==Dir::imax(), 
       j = bc.dir==Dir::jmin() || bc.dir==Dir::jmax(),
       k = bc.dir==Dir::kmin() || bc.dir==Dir::kmax();
  bool b = bc.dir==Dir::ibody();                              

  if( j || k || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.ir.first() < cxg.first()+boil::BW-1 && bc.ir.last() < cxg.first()+boil::BW-1 ||
        bc.ir.first() > cxg.last()+boil::BW-1 && bc.ir.last() > cxg.last()+boil::BW-1 ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(cxg.contains(bc.ir.first()-boil::BW+1)) bc.ir.first(bc.ir.first()-cxg.first()+1);
      else                          bc.ir.first(1);
      if(cxg.contains(bc.ir.last()-boil::BW+1)) bc.ir.last(bc.ir.last()-cxg.first()+1);
      else                         bc.ir.last(dom->ni()-boil::BW-1);
    }
  } /* j || k || b */

  if( i || k || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.jr.first() < cyg.first()+boil::BW-1 && bc.jr.last() < cyg.first()+boil::BW-1 ||
        bc.jr.first() > cyg.last()+boil::BW-1 && bc.jr.last() > cyg.last()+boil::BW-1 ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(cyg.contains(bc.jr.first()-boil::BW+1)) bc.jr.first(bc.jr.first()-cyg.first()+1);
      else                          bc.jr.first(1);
      if(cyg.contains(bc.jr.last()-boil::BW+1)) bc.jr.last(bc.jr.last()-cyg.first()+1);
      else                          bc.jr.last(dom->nj()-boil::BW-1);
    }
  } /* i || k || b */

  if( i || j || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.kr.first() < czg.first()+boil::BW-1 && bc.kr.last() < czg.first()+boil::BW-1 ||
        bc.kr.first() > czg.last()+boil::BW-1 && bc.kr.last() > czg.last()+boil::BW-1 ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(czg.contains(bc.kr.first()-boil::BW+1)) bc.kr.first(bc.kr.first()-czg.first()+1);
      else                          bc.kr.first(1);
      if(czg.contains(bc.kr.last()-boil::BW+1)) bc.kr.last(bc.kr.last()-czg.first()+1);
      else                          bc.kr.last(dom->nk()-boil::BW-1);
    }
  } /* i || j || b */
#else
  /*-------------------------------------------------+
  |  2. set range for directions explicitly defined  | 
  +-------------------------------------------------*/

  /* this is correct: in the aligned direction,
     the idx should be one just outside the domain in the given dir */
  if(bc.dir == Dir::imin()) {bc.ir.first(boil::BW-1);
                             bc.ir.last (boil::BW-1);}
  if(bc.dir == Dir::imax()) {bc.ir.first(dom->ni()-boil::BW);
                             bc.ir.last (dom->ni()-boil::BW);}

  if(bc.dir == Dir::jmin()) {bc.jr.first(boil::BW-1);
                             bc.jr.last (boil::BW-1);}
  if(bc.dir == Dir::jmax()) {bc.jr.first(dom->nj()-boil::BW);
                             bc.jr.last (dom->nj()-boil::BW);}

  if(bc.dir == Dir::kmin()) {bc.kr.first(boil::BW-1);
                             bc.kr.last (boil::BW-1);}
  if(bc.dir == Dir::kmax()) {bc.kr.first(dom->nk()-boil::BW);
                             bc.kr.last (dom->nk()-boil::BW);}

  /*--------------------------------------+
  |  3. set the range for the directions  |
  |         not explicitly defined        |
  +--------------------------------------*/
  Range<int> cxg = dom->cxg(); 
  Range<int> cyg = dom->cyg(); 
  Range<int> czg = dom->czg(); 

  /* c*g contains the indices on the given processor
     on the global domain, w/o halo cells
     e.g. (33,64) out of (1,64)                     */

  if( bc.ir.first() == 0 && bc.ir.last() == -1 ) {
    bc.ir.first(cxg.first()); // as sx()
    bc.ir.last (cxg.last() ); // as ex()
  }
  if( bc.jr.first() == 0 && bc.jr.last() == -1 ) {
    bc.jr.first(cyg.first()); // as sy()
    bc.jr.last (cyg.last() ); // as ey()
  }
  if( bc.kr.first() == 0 && bc.kr.last() == -1 ) {
    bc.kr.first(czg.first()); // as sz()
    bc.kr.last (czg.last() ); // as ez()
  }

  /* now the non-aligned indices correspond exactly to the
     c*g, provided that they were not user-specified 
     and they have to be corrected = transformed to local
     indices, accounting for halo cells */

  /*---------------------------------------------------+
  |  4. correct the range for directions defined by 3  |
  +---------------------------------------------------*/
  bool i = bc.dir==Dir::imin() || bc.dir==Dir::imax(), 
       j = bc.dir==Dir::jmin() || bc.dir==Dir::jmax(),
       k = bc.dir==Dir::kmin() || bc.dir==Dir::kmax();
  bool b = bc.dir==Dir::ibody();                              

  if( j || k || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.ir.first() < cxg.first() && bc.ir.last() < cxg.first() ||
        bc.ir.first() > cxg.last() && bc.ir.last() > cxg.last() ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(cxg.contains(bc.ir.first())) bc.ir.first(bc.ir.first()-cxg.first()+boil::BW);
      else                          bc.ir.first(boil::BW);
      if(cxg.contains(bc.ir.last())) bc.ir.last(bc.ir.last()-cxg.first()+boil::BW);
      else                         bc.ir.last(dom->ni()-boil::BW-1);
    }
  } /* j || k || b */

  if( i || k || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.jr.first() < cyg.first() && bc.jr.last() < cyg.first() ||
        bc.jr.first() > cyg.last() && bc.jr.last() > cyg.last() ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(cyg.contains(bc.jr.first())) bc.jr.first(bc.jr.first()-cyg.first()+boil::BW);
      else                          bc.jr.first(boil::BW);
      if(cyg.contains(bc.jr.last())) bc.jr.last(bc.jr.last()-cyg.first()+boil::BW);
      else                          bc.jr.last(dom->nj()-boil::BW-1);
    }
  } /* i || k || b */

  if( i || j || b ) {

    /* if both outside the range - skip this b.c. */
    if( bc.kr.first() < czg.first() && bc.kr.last() < czg.first() ||
        bc.kr.first() > czg.last() && bc.kr.last() > czg.last() ) {
      bc.ir.first( 0); bc.jr.first( 0); bc.kr.first( 0);
      bc.ir.last (-1); bc.jr.last (-1); bc.kr.last (-1);
    }
    else {
      if(czg.contains(bc.kr.first())) bc.kr.first(bc.kr.first()-czg.first()+boil::BW);
      else                          bc.kr.first(boil::BW);
      if(czg.contains(bc.kr.last())) bc.kr.last(bc.kr.last()-czg.first()+boil::BW);
      else                          bc.kr.last(dom->nk()-boil::BW-1);
    }
  } /* i || j || b */
#endif

  //boil::oout<<"Bnd after 4 "<<bc.dir<<" "<<bc.ir.first()<<" "<<bc.ir.last()<<" | "<<bc.kr.first()<<" "<<bc.kr.last()<<" || "<<cxg.first()<<" "<<cxg.last()<<" | "<<czg.first()<<" "<<czg.last()<<boil::endl;

  section.push_back(bc);

  /*-------------------------------------------------------------------+
  |  Note: loop range for vector used for Momentum object is reset in: |
  |  staggered_set_ranges.cpp                                          | 
  +-------------------------------------------------------------------*/
}

/******************************************************************************/
BndType & BndCnd::type(const int b) {
  return section.at(b).typ;
}

/******************************************************************************/
char * BndCnd::formula(const int b, const Comp m) const {
  return section.at(b).val[~m].c;
}

/******************************************************************************/
real BndCnd::value(const int b, const Comp m) const {
  return section.at(b).val[~m].r;
}

/******************************************************************************/
real BndCnd::value(const int b, const real x, const real y, const real z, 
                   const Comp m) const {

  /* if formula is not defined */
  if(!section.at(b).val[~m].c)
    return section.at(b).val[~m].r;

  /* if formula is defined */
  else {
    Formula f;

    std::stringstream xs, ys, zs;
    xs << "x=" << x; f.evaluate(xs);
    ys << "y=" << y; f.evaluate(ys);
    zs << "z=" << z; f.evaluate(zs);

    return f.evaluate(section.at(b).val[~m].c);
  }
}

/******************************************************************************/
int BndCnd::count() const {
  return section.size();
}

/******************************************************************************/
const Dir & BndCnd::direction(const int b) const {
  return section.at(b).dir;
}

/******************************************************************************/
bool BndCnd::type(const Dir & dir, const BndType & bt) const {

  for(int b=0; b<count(); b++) 
    if( section.at(b).typ == bt &&              
        section.at(b).dir == dir ) return true;

  return false;
}

/******************************************************************************/
int BndCnd::index(const Dir & dir, const BndType & bt) const {

  for(int b=0; b<count(); b++){
    if( section.at(b).typ == bt &&
        section.at(b).dir == dir ) return b;
  }
  return -1;
}

/******************************************************************************/
bool BndCnd::type_here(const Dir & dir, const BndType & bt) const {

  if( dir == Dir::imin() && dom->coord(Comp::i()) != 0 )             
    return false;
  if( dir == Dir::imax() && dom->coord(Comp::i()) != dom->dim(Comp::i())-1 ) 
    return false;
  if( dir == Dir::jmin() && dom->coord(Comp::j()) != 0 )             
    return false;
  if( dir == Dir::jmax() && dom->coord(Comp::j()) != dom->dim(Comp::j())-1 ) 
    return false;
  if( dir == Dir::kmin() && dom->coord(Comp::k()) != 0 )             
    return false;
  if( dir == Dir::kmax() && dom->coord(Comp::k()) != dom->dim(Comp::k())-1 ) 
    return false;

  for(int b=0; b<count(); b++) 
    if( section.at(b).typ == bt &&              
        section.at(b).dir == dir ) return true;

  return false;
}

/******************************************************************************/
bool BndCnd::type_decomp(const int b) const {
  
  Dir dir = section.at(b).dir;
  
  if( dir == Dir::imin() && dom->coord(Comp::i()) == 0 )
    return false;
  if( dir == Dir::imax() && dom->coord(Comp::i()) == dom->dim(Comp::i())-1 )
    return false;
  if( dir == Dir::jmin() && dom->coord(Comp::j()) == 0 )
    return false;
  if( dir == Dir::jmax() && dom->coord(Comp::j()) == dom->dim(Comp::j())-1 )
    return false;
  if( dir == Dir::kmin() && dom->coord(Comp::k()) == 0 )
    return false;
  if( dir == Dir::kmax() && dom->coord(Comp::k()) == dom->dim(Comp::k())-1 )
    return false;
  
  return true;
}

/******************************************************************************/
bool BndCnd::type_decomp(const Dir dir) const {

  if( dir == Dir::imin() && dom->coord(Comp::i()) == 0 )
    return false;
  if( dir == Dir::imax() && dom->coord(Comp::i()) == dom->dim(Comp::i())-1 )
    return false;
  if( dir == Dir::jmin() && dom->coord(Comp::j()) == 0 )
    return false;
  if( dir == Dir::jmax() && dom->coord(Comp::j()) == dom->dim(Comp::j())-1 )
    return false;
  if( dir == Dir::kmin() && dom->coord(Comp::k()) == 0 )
    return false;
  if( dir == Dir::kmax() && dom->coord(Comp::k()) == dom->dim(Comp::k())-1 )
    return false;

  return true;
}

/******************************************************************************/
bool BndCnd::exists(const BndType & bt) const {

  int here=0;

  /* look for bt among boundary conditions */
  for(int b=0; b<count(); b++) {
    if( section.at(b).typ == bt ) {here=1; break;}
  }

  /* check if it is found on other processor(s) */
  boil::cart.sum_int(&here);

  if(here) 
    return true;

  return false;
}

/******************************************************************************/
bool BndCnd::exists(const int b) const {
  return section.at(b).ir.exists() &&
         section.at(b).jr.exists() &&
         section.at(b).kr.exists();
}

/***************************************************************************//**
*  Prints the name of the current booundary condition. 
*  Needed only for information or debugging.
*******************************************************************************/
std::ostream & operator << (std::ostream &ost, const BndCnd & bc) {

  ost << std::endl;
  ost << "type: " << bc.typ << std::endl;
  ost << "dir:  " << bc.dir << std::endl;
  ost << "i: "    << bc.ir  << std::endl;
  ost << "j: "    << bc.jr  << std::endl;
  ost << "k: "    << bc.kr;

  return ost;
}

/***************************************************************************//**
*  Prints the name of the current booundary type. 
*  Needed only for information or debugging.
*******************************************************************************/
std::ostream & operator << (std::ostream &ost, const BndType & bt) {

  switch(bt.val) {
    case( 0): ost << "undefined ";  break;
    case( 1): ost << "dirichlet ";  break;
    case( 2): ost << "neumann ";    break;
    case( 3): ost << "periodic ";   break;
    case( 4): ost << "inlet ";      break;
    case( 5): ost << "outlet ";     break;
    case( 6): ost << "wall ";       break;
    case( 7): ost << "symmetry ";   break;
    case( 8): ost << "insert ";     break;
    case( 9): ost << "convective "; break;
  }
    
  return ost;
}
