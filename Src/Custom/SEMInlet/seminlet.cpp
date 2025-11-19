// =========================== seminlet.cpp ===========================
#include "seminlet.h"
#include <fstream>
#include <sstream>

// --- ctor: nothing configured yet
SEMInlet::SEMInlet(uint64_t seed)
: LY_(0), LZ_(0), sig_rep_(0), Uc_(0),
  B_{0,0,0,0}, rng_(seed), ready_(false), ifunc_(1) {} 
  // ifunc_(1) means tent is used as a default

// --- configure geometry & parameters (no eddies yet)
void SEMInlet::setup(real Y_min, real Y_max, real Z_min, real Z_max,
                     real sigma_rep, real Uc, real alp) {
  if (Y_max-Y_min<=0 || Z_max-Z_min<=0 || sigma_rep<=0 || Uc<=0)
    throw std::runtime_error("SEMInlet::setup: non-positive parameter");
  LY_ = Y_max-Y_min; LZ_ = Z_max-Z_min; sig_rep_ = sigma_rep; Uc_ = Uc;
  B_ = {Y_min, Y_max, Z_min, Z_max};
  alp_ = alp;
  real LX = 2.0 * sigma_rep;                        // -sigma_rep < X < sigma_rep
  n_eddy_ = alp * (LX*LY_*LZ_) / (pow(sigma_rep,3));  // Domain volume: LX*LY*LZ
  ready_ = true; // configured (no eddies yet)
}

// --- optionally change bounds after setup
void SEMInlet::setBounds(real y_min, real y_max, real z_min, real z_max) {
  B_ = {y_min, y_max, z_min, z_max};
}

// --- initialize/resize eddies
void SEMInlet::initEddies() {
  requireReady();
  E_.assign(n_eddy_, Eddy{});
  for (auto &e : E_) { randomizeEddy(e, /*placeAnywhere=*/true); }
}

// --- reseed RNG
void SEMInlet::setSeed(uint64_t seed) { rng_.seed(seed); }

// --- advance dt & recycle at +sig_rep
void SEMInlet::advance(real dt) {
  requireReady();
  for (auto &e : E_) {
    e.x += Uc_ * dt;
    if (e.x > +sig_rep_) recycle(e);
  }
}

// --- accumulate raw SEM signals s_i at (y,z)
void SEMInlet::accumulate(real y, real z, real s_raw[3]) const {
  requireReady();
  s_raw[0] = s_raw[1] = s_raw[2] = 0.0;
  //const real inv_sqrtN = 1.0 / std::sqrt(real(E_.size()));
  // sqrt(Vb/(n_eddy*sig^3))=sqrt(1/alp) Eq.(4) Schau (2022) C&F 105671
  // https://doi.org/10.1016/j.compfluid.2022.105671
  const real inv_sqrtN = 1.0 / std::sqrt(alp_);
  for (const auto& e : E_) {
    const real dx = -e.x;       // inlet plane x=0
    const real dy =  y - e.y0;
    const real dz =  z - e.z0;
    real kern;
    if (ifunc_==0) {
      kern = gauss_iso(dx, dy, dz);
    } else if (ifunc_==1) {
      kern = tent(dx, dy, dz);
    } else {
      boil::oout<<"#Error!!! seminlet.cpp exiting!\n";
      exit(0);
    }
    s_raw[0] += real(e.sgn[0]) * kern;
    s_raw[1] += real(e.sgn[1]) * kern;
    s_raw[2] += real(e.sgn[2]) * kern;
  }
  s_raw[0] *= inv_sqrtN;
  s_raw[1] *= inv_sqrtN;
  s_raw[2] *= inv_sqrtN;
}

// --- private helpers
void SEMInlet::randomizeEddy(Eddy& e, bool placeAnywhere) {
  // placeAnywhere=true → x ∈ [-sig_rep, +sig_rep] at start; recycle → x = -sig_rep
  e.x  = placeAnywhere ? urand(rng_, -sig_rep_, +sig_rep_) : -sig_rep_;
  e.y0 = urand(rng_, B_.y_min, B_.y_max);
  e.z0 = urand(rng_, B_.z_min, B_.z_max);
  for (int i=0;i<3;++i) e.sgn[i] = rand_sign(rng_);
}

void SEMInlet::recycle(Eddy& e) {
  randomizeEddy(e, /*placeAnywhere=*/false);
}

// --- plot for visualization
void SEMInlet::plot(){
  if (boil::cart.iam()==0) {
    std::ofstream fout("SEMInlet.dat");
    // Tecplot ASCII header
    fout << "TITLE     = \"Scattered 3D Points\"\n";
    fout << "VARIABLES = \"X\" \"Y\" \"Z\"\n";
    fout << "ZONE T=\"Scatter\", I=" << E_.size() << ", F=POINT\n";
    // Write point data
    for (auto &e : E_) {
      fout << e.x << " " << e.y0 << " " << e.z0 << "\n";
    }
    fout.close();
    boil::oout << "# Plotting: sem.dat\n";
  }
}

// --- save bck files
void SEMInlet::save(const char * nm, const int it) {
  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());
  /* open a file */
  std::ofstream out(name.c_str(), std::ios::binary);
  /* save necessary variables */
    // RNG state (for random number)
    auto state = rng_; // copy
    std::ostringstream oss;
    oss << state; // serialize engine to text
    std::string rng_str = oss.str();
    size_t rng_size = rng_str.size();
    out.write(reinterpret_cast<char*>(&rng_size), sizeof(rng_size));
    out.write(rng_str.data(), rng_size);
    // eddy vector
    size_t N = E_.size();
    out.write(reinterpret_cast<char*>(&N), sizeof(N));
    out.write(reinterpret_cast<const char*>(E_.data()), N*sizeof(Eddy));
  /* close a file */
  out.close();
}

// --- load bck files
void SEMInlet::load(const char * nm, const int it) {
  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());
  /* open a file */
  std::ifstream in(name.c_str(), std::ios::binary);
  /* load variables */
    if(!in) throw std::runtime_error("SEMInlet::load: cannot open file");
    // RNG state
    size_t rng_size;
    in.read(reinterpret_cast<char*>(&rng_size), sizeof(rng_size));
    std::string rng_str(rng_size, '\0');
    in.read(&rng_str[0], rng_size);
    std::istringstream iss(rng_str);
    iss >> rng_;
    // eddy vector
    size_t N;
    in.read(reinterpret_cast<char*>(&N), sizeof(N));
    E_.resize(N);
    in.read(reinterpret_cast<char*>(E_.data()), N*sizeof(Eddy));
  /* close a file */
  in.close();
}

void SEMInlet::rm(const char * nm, const int it) {
  /* file name */
  std::string name = name_file(nm, ".bck", it, boil::cart.iam());
  /* open a file */
  remove(name.c_str());
}

