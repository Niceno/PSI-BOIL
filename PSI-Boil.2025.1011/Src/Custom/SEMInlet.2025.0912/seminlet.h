#ifndef SEMINLET_H
#define SEMINLET_H
// =========================== seminlet.h ===========================
// Original SEM (Jarrin 2006)  symmetric slab, NO EWMA
// Axes: x(streamwise), y(spanwise), z(wall-normal)

#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include "../../Global/global_name_file.h"
#include "../../Parallel/communicator.h"
#include "../../Parallel/Out/print.h"

class SEMInlet {
public:
  struct BoundsYZ { real y_min, y_max, z_min, z_max; };
  struct Eddy {
    real x, y0, z0;   // eddy center
    int8_t sgn[3];      // independent signs for {sx, sy, sz}
  };
  //struct Lchol { real L11, L22, L31, L33; };

  // --- ctor (lightweight; call setup() then initEddies())
  explicit SEMInlet(uint64_t seed = 1234567ULL);

  // --- configure geometry & params (no eddies yet)
  void setup(real LY, real LZ, real sigma_rep, real Uc);

  // --- optionally override bounds (after setup)
  void setBounds(real y_min, real y_max, real z_min, real z_max);

  // --- initialize/resize eddies
  void initEddies(int N);

  // --- reseed RNG
  void setSeed(uint64_t seed);

  // --- advance eddies dt and recycle at +sigma_rep
  void advance(real dt);

  // --- accumulate raw SEM signals s_i at inlet node (y,z)
  //     using caller-provided local sigma(z)
  void accumulate(real y, real z, real sigma_local, real s_raw[3]) const;

  // --- save ,load and rm
  void save(const char *, const int);
  void load(const char *, const int);
  void rm  (const char *, const int);

  // --- accessors
  int    size()     const { return int(E_.size()); }
  real sigmaRep() const { return sig_rep_; }
  real Uc()       const { return Uc_; }
  BoundsYZ bounds() const { return B_; }

  // Jarrin Gaussian prefactor: sqrt(8/pi^(3/2))
  static constexpr real Kpref = 1.5163266492815837;
 
  void set_func(const std::string& name) {
    if (name=="gauss") {
      ifunc_=0;
      boil::oout<<"seminlet::use gauss\n";
    } else if (name=="tent") {
      ifunc_=1;
      boil::oout<<"seminlet::use tent\n";
    } else {
      boil::oout<<"seminlet::set_func Error!!! gaus or tent. exiting! \n";
      exit(0);
    }
  }

private:
  // helpers (inline/private)
  static inline int8_t rand_sign(std::mt19937_64& rng) {
    std::uniform_int_distribution<int> d(0,1);
    return d(rng) ? int8_t(1) : int8_t(-1);
  }
  static inline real urand(std::mt19937_64& rng, real a, real b) {
    std::uniform_real_distribution<real> d(a,b);
    return d(rng);
  }
  static inline real gauss_iso(real dx, real dy, real dz, real sigma) {
    const real inv2s2 = 1.0 / (2.0 * sigma * sigma);
    return Kpref * std::exp( - (dx*dx + dy*dy + dz*dz) * inv2s2 );
  }
  static inline real tent(real dx, real dy, real dz, real sigma) {
    return  sqrt(3.0/2.0)*std::max(0.0,(1.0-std::abs(dx/sigma)))
           *sqrt(3.0/2.0)*std::max(0.0,(1.0-std::abs(dy/sigma)))
           *sqrt(3.0/2.0)*std::max(0.0,(1.0-std::abs(dz/sigma)));
  }

  void randomizeEddy(Eddy& e, bool placeAnywhere);
  void recycle(Eddy& e);
  void requireReady() const {
    if (!ready_) throw std::runtime_error("SEMInlet: call setup(...) first");
  }

  // data
  real LY_, LZ_;
  real sig_rep_;
  real Uc_;
  BoundsYZ B_;
  std::mt19937_64 rng_;
  std::vector<Eddy> E_;
  bool ready_;
  int ifunc_;
};

#endif
