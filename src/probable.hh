#pragma once

#include <cstdint>
#include <cxxabi.h>
#include <ostream>
#include <typeinfo>
#include <vector>

#include <units.hh>
#include <vec3.hh>

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

class Dumper {
  virtual void dump(int particle, int step, const Vec3 &r, const Vec3 &p) = 0;
};

void simulate(const std::vector<Scattering *> mechanisms,
              double temperature,
              const Vec3 &electric_field,
              const Vec3 &magnetic_field,
              double time_step,
              double all_time,
              size_t ensemble_size,
              Dumper &dumper);

// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable
