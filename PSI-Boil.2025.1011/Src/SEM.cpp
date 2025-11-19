// =========================== SEM.hpp ===========================
// Original SEM (Jarrin 2006)  symmetric slab, NO EWMA normalization
// Axes: x(streamwise), y(spanwise), z(wall-normal)

#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <cstdint>

// --------- Types ----------
struct BoundsYZ {
  double y_min, y_max;  // spanwise range
  double z_min, z_max;  // wall-normal range
};

struct Eddy {
  double x;    // streamwise eddy center in [-sig_rep, +sig_rep]
  double y0;   // spanwise eddy center
  double z0;   // wall-normal eddy center
  int8_t sgn[3]; // independent signs for {sx, sy, sz} contributions
};

struct Lchol {
  // Cholesky entries for R = [[Rxx,0,Rxz],[0,Ryy,0],[Rxz,0,Rzz]]
  double L11, L22, L31, L33;
};

// --------- RNG helpers ----------
inline int8_t rand_sign(std::mt19937_64& rng){
  std::uniform_int_distribution<int> d(0,1);
  return d(rng) ? int8_t(1) : int8_t(-1);
}
inline double urand(std::mt19937_64& rng, double a, double b){
  std::uniform_real_distribution<double> d(a,b);
  return d(rng);
}

// --------- Gaussian kernel (isotropic) ---------
// Jarrin Eq.(8) with SEM normalization: sqrt(8 / pi^(3/2))
constexpr double Kpref = 1.5163266492815837; // = sqrt(8.0 / pow(pi,1.5))
inline double gauss_iso(double dx, double dy, double dz, double sigma_local){
  const double inv2s2 = 1.0 / (2.0 * sigma_local * sigma_local);
  return Kpref * std::exp( - (dx*dx + dy*dy + dz*dz) * inv2s2 );
}

// --------- Initialize eddies in symmetric slab ---------
inline void sem_init_eddies(std::vector<Eddy>& E, int N,
                            const BoundsYZ& B, double sig_rep,
                            std::mt19937_64& rng)
{
  E.resize(N);
  for (int n=0; n<N; ++n) {
    E[n].x  = urand(rng, -sig_rep, +sig_rep);
    E[n].y0 = urand(rng, B.y_min, B.y_max);
    E[n].z0 = urand(rng, B.z_min, B.z_max);
    for (int i=0;i<3;++i) E[n].sgn[i] = rand_sign(rng);
  }
}

// --------- Convect and recycle (only at +x boundary) ---------
inline void sem_advect_recycle(std::vector<Eddy>& E,
                               double Uc, double dt,
                               const BoundsYZ& B,
                               double sig_rep,
                               std::mt19937_64& rng)
{
  for (auto& e : E) {
    e.x += Uc * dt;
    if (e.x > +sig_rep) {            // crossed the +x boundary → recycle upstream
      e.x  = -sig_rep;
      e.y0 = urand(rng, B.y_min, B.y_max);
      e.z0 = urand(rng, B.z_min, B.z_max);
      for (int i=0;i<3;++i) e.sgn[i] = rand_sign(rng);
    }
  }
}

// --------- Accumulate raw SEM signals s_i at (y,z) ---------
// NOTE: we use sigma_local = sigma(z_target) (no per-eddy sigma).
inline void sem_accumulate_raw(const std::vector<Eddy>& E,
                               double y, double z,
                               double sigma_local,
                               double s_raw[3])
{
  s_raw[0] = s_raw[1] = s_raw[2] = 0.0;
  const double inv_sqrtN = 1.0 / std::sqrt(double(E.size()));

  for (const auto& e : E) {
    const double dx = -e.x;            // inlet plane at x=0
    const double dy =  y - e.y0;
    const double dz =  z - e.z0;
    const double kern = gauss_iso(dx, dy, dz, sigma_local);
    s_raw[0] += double(e.sgn[0]) * kern;
    s_raw[1] += double(e.sgn[1]) * kern;
    s_raw[2] += double(e.sgn[2]) * kern;
  }
  s_raw[0] *= inv_sqrtN;
  s_raw[1] *= inv_sqrtN;
  s_raw[2] *= inv_sqrtN;
}

// --------- Map to target stresses with Cholesky ---------
inline void sem_apply_cholesky(const Lchol& L,
                               const double s_raw[3],
                               double& up, double& vp, double& wp)
{
  // u' = L11*sx
  up = L.L11 * s_raw[0];
  // v' = L22*sy
  vp = L.L22 * s_raw[1];
  // w' = L31*sx + L33*sz
  wp = L.L31 * s_raw[0] + L.L33 * s_raw[2];
}

