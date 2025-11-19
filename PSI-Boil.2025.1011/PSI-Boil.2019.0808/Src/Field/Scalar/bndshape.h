/***************************************************************************//**
*  \brief Ravioli class for boundary flagging for faster if-condition evals.
*******************************************************************************/
#ifndef BNDFLAG_H
#define BNDFLAG_H

/* this should be templated for all Field objects! */

///////////////
//           //
//  BndFlag  //
//           //
///////////////
class BndFlag {
  public:
    BndFlag(const Scalar & sca) :
      iminp(false), imaxp(false), /* true for periodic */ 
      jminp(false), jmaxp(false),
      kminp(false), kmaxp(false),
      iminc(true), imaxc(true), /* true for cut-stencil */ 
      jminc(true), jmaxc(true),
      kminc(true), kmaxc(true),
      iminw(false), imaxw(false), /* true for wall */ 
      jminw(false), jmaxw(false),
      kminw(false), kmaxw(false),
      ifull(true), jfull(true), kfull(true) /* true for not a dummy direction */
  {
    // imin
    Dir d = Dir::imin();
    if (sca.bc().type_decomp(d)) {
      iminp=true;
      iminc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        iminp=true;
        iminc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        iminw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        iminp=true;
        iminc=false;
        ifull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) iminc=false;
    }
    // imax
    d = Dir::imax();
    if (sca.bc().type_decomp(d)) {
      imaxp=true;
      imaxc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        imaxp=true;
        imaxc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        imaxw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        imaxp=true;
        imaxc=false;
        ifull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) imaxc=false;
    }
    // jmin
    d = Dir::jmin();
    if (sca.bc().type_decomp(d)) {
      jminp=true;
      jminc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        jminp=true;
        jminc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        jminw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        jminp=true;
        jminc=false;
        jfull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) jminc=false;
    }
    // jmax
    d = Dir::jmax();
    if (sca.bc().type_decomp(d)) {
      jmaxp=true;
      jmaxc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        jmaxp=true;
        jmaxc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        jmaxw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        jmaxp=true;
        jmaxc=false;
        jfull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) jmaxc=false;
    }
    // kmin
    d = Dir::kmin();
    if (sca.bc().type_decomp(d)) {
      kminp=true;
      kminc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        kminp=true;
        kminc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        kminw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        kminp=true;
        kminc=false;
        kfull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) kminc=false;
    }
    // kmax
    d = Dir::kmax();
    if (sca.bc().type_decomp(d)) {
      kmaxp=true;
      kmaxc=false;
    } else {
      if (sca.bc().type(d,BndType::periodic())) {
        kmaxp=true;
        kmaxc=false;
      } else if (sca.bc().type(d,BndType::wall())) {
        kmaxw=true;
      } else if (sca.bc().type(d,BndType::pseudo())) {
        kmaxp=true;
        kmaxc=false;
        kfull = false;
      }
      if (sca.domain()->bnd_symmetry(d)) kmaxc=false;
    }
 
    dim = int(ifull)+int(jfull)+int(kfull);

  } /* constructor */

  bool iminp, imaxp, /* true for periodic */
       jminp, jmaxp,
       kminp, kmaxp,
       iminc, imaxc, /* true for cut-stencil */
       jminc, jmaxc,
       kminc, kmaxc,
       iminw, imaxw, /* true for wall */
       jminw, jmaxw,
       kminw, kmaxw,
       ifull, jfull, kfull;
  int dim; /* true dimension of the problem */
};
#endif
