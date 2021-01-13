#pragma once

#include <cstdint>
#include <cxxabi.h>
#include <ostream>
#include <typeinfo>
#include <vector>

#include <vec3.hh>
#include <units.hh>

namespace probable {

struct Band {
  bool occupied;

  virtual double energy(const Vec3 &p) const = 0;
  virtual Vec3 velocity(const Vec3 &p) const = 0;

  virtual bool try_boltzmann_sample_momentum(double temperature, Vec3 &momentum) const = 0;

  virtual double acoustic_scattering_integral(double energy) const = 0;
  virtual double optical_scattering_integral(double energy) const = 0;
  virtual Vec3 acoustic_scatter(double energy, double random) const = 0;
  virtual Vec3 optical_scatter(double energy, double random) const = 0;

  Band(bool occupied) : occupied(occupied) {}
  virtual ~Band() {}
};

struct Particle {
  Vec3 r;
  Vec3 p;
  Band &band;
};

struct Material {
  double acoustic_deformation_potential;
  double nonpolar_optical_deformation_potential;
  double sound_velocity;
  double density;

  std::vector<Band *> bands;
  Particle create_particle(double temperature) const;
  ~Material() {
    for (auto b : bands) {
      delete b;
    }
  }
};

struct ParabolicBand : public Band {
  double mass;
  double energy(const Vec3 &p) const { return p.dot(p) / (2 * mass); }
  Vec3 velocity(const Vec3 &p) const { return p / mass; }
  bool try_boltzmann_sample_momentum(double temperature, Vec3 &momentum) const;
  double acoustic_scattering_integral(double energy) const { return 0; }
  double optical_scattering_integral(double energy) const { return 0; }
  Vec3 acoustic_scatter(double energy, double random) const { return Vec3(); }
  Vec3 optical_scatter(double energy, double random) const { return Vec3(); }
  ParabolicBand(bool occupied, double mass) : Band(occupied), mass(mass) {}
};

struct Scattering {
  const Material &material;
  const Band &band;
  const double energy;
  Scattering(const Material &m, const Band &b, double e) : material(m), band(b), energy(e) {}
  virtual double rate(const Vec3 &p) const = 0;
  virtual Vec3 scatter(const Vec3 &p) const = 0;
  virtual ~Scattering() {}
};

struct AcousticScattering : public Scattering {
  AcousticScattering(const Material &m, const Band &b) : Scattering(m, b, 0) {}
  virtual double rate(const Vec3 &p) const {
    return pow(material.acoustic_deformation_potential, 2) * consts::kB * temperature *
           (math::pi * pow(consts::hbar, 4) * material.density * pow(material.sound_velocity, 2)) *
           band.acoustic_scattering_integral(band.energy(p));
  }
  virtual Vec3 scatter(const Vec3 &p, double random) const {
    return band.acoustic_scatter(band.energy(p), random);
  };
};

struct OpticalScattering : public Scattering {
  OpticalScattering(const Material &m, const Band &b, double energy) : Scattering(m, b, energy) {}
  virtual double rate(const Vec3 &p) const {
    return band.optical_scattering_integral(band.energy(p) - energy);
  }
  virtual Vec3 scatter(const Vec3 &p, double random) const {
    return band.optical_scatter(band.energy(p) - energy, random);
  };
};

inline std::ostream &operator<<(std::ostream &s, const Scattering &sc) {
  int status;
  char *realname = abi::__cxa_demangle(typeid(sc).name(), 0, 0, &status);
  s << realname << " " << sc.energy;
  free(realname);
  return s;
}

enum DumpFlags {
  // contents
  none = 0,
  number = 1,
  time = number << 1,
  momentum = time << 1,
  energy = momentum << 1,
  velocity = energy << 1,
  scattering = velocity << 1,
  all = number | time | momentum | energy | velocity | scattering,
  // frequency
  // without this flag it will dump on every step
  on_scatterings = scattering << 1,
};

struct Results {
  size_t size;
  DumpFlags flags;
  Vec3 average_velocity;
  std::vector<uint32_t> scattering_count;
  std::vector<uint32_t> ns;
  std::vector<double> ts;
  std::vector<Vec3> momentums;
  std::vector<Vec3> velocities;
  std::vector<double> energies;
  std::vector<uint32_t> scatterings;
  Results() {}
  Results(size_t cap, DumpFlags flags = DumpFlags::none)
      : size(0), flags(flags), ns(), ts(), momentums(), velocities(), energies(), scatterings() {
    ns.reserve(cap);
    ts.reserve(cap);
    momentums.reserve(cap);
    velocities.reserve(cap);
    energies.reserve(cap);
    scatterings.reserve(cap);
  }
  void append(uint32_t n, double t, const Vec3 &p, const Vec3 &v, double e, size_t s);
  friend std::ostream &operator<<(std::ostream &s, const Results &r);
};

std::vector<Results> simulate(const std::vector<Scattering *> mechanisms,
                              double temperature,
                              const Vec3 &electric_field,
                              const Vec3 &magnetic_field,
                              double time_step,
                              double all_time,
                              size_t ensemble_size,
                              DumpFlags flags = DumpFlags::all);

// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable
