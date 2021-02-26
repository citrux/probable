#pragma once

#include <cstdint>
#include <cxxabi.h>
#include <functional>
#include <ostream>
#include <typeinfo>
#include <vector>

#include <units.hh>
#include <vec3.hh>

namespace probable {

struct Band {
  bool occupied;

  /// dispersion relation for band
  virtual double energy(const Vec3 &p) const = 0;
  /// d/dp energy
  virtual Vec3 velocity(const Vec3 &p) const = 0;

  /// try to create new particle, true if created, false if not
  virtual bool try_boltzmann_sample_momentum(double temperature, Vec3 &momentum) const = 0;

  /// \int_{\epsilon(p)=energy} \frac{d\sigma}{|\nabla\varepsilon(\vec{p})|}
  virtual double delta_integral(double energy) const = 0;

  /// deprecated
  /// optimizations for acoustic scattering
  /// scattering integral (see paper)
  virtual double acoustic_scattering_integral(double energy) const = 0;
  /// new momentum after scattering
  virtual Vec3 acoustic_scatter(double energy, std::function<double()> random) const = 0;

  /// optimizations for optical scattering
  /// scattering integral (see paper)
  /// \int_{\epsilon(p)=energy} \frac{(\vec{p} - \vec{p}_0)^2}{((\vec{p} - \vec{p}_0)^2 + \hbar^2 q_0^2)^2}\frac{d\sigma}{|\nabla\varepsilon(\vec{p})|}
  virtual double optical_scattering_integral(double energy, const Vec3 &momentum, double phonon_energy) const = 0;
  /// new momentum after scattering
  virtual Vec3 optical_scatter(double energy, std::function<double()> random) const = 0;

  Band(bool occupied) : occupied(occupied) {}
  virtual ~Band() {}
};

struct Particle {
  Vec3 r;
  Vec3 p;
  double charge;
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
  double delta_integral(double energy) const { return sqrt(2 * mass * energy); }
  double acoustic_scattering_integral(double energy) const { return delta_integral(energy); }
  double optical_scattering_integral(double energy, const Vec3 &momentum, double phonon_energy) const {
    return mass / momentum.length() * asinh(energy / phonon_energy);
  }
  Vec3 acoustic_scatter(double energy, std::function<double()> random) const {
    double p = sqrt(2 * mass * energy);
    double r = random();
    double cos_theta = 1 - 2 * r;
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = 2 * math::pi * random();
    return p * Vec3{sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
  }
  Vec3 optical_scatter(double energy, std::function<double()> random) const { return acoustic_scatter(energy, random); }
  ParabolicBand(bool occupied, double mass) : Band(occupied), mass(mass) {}
};

struct Scattering {
  const Material &material;
  const Band &band;
  const double energy;

  Scattering(const Material &m, const Band &b, double e) : material(m), band(b), energy(e) {}
  virtual double rate(const Vec3 &p) const = 0;
  virtual Vec3 scatter(const Vec3 &p, std::function<double()> random) const = 0;
  virtual ~Scattering() {}
};

struct AcousticScattering : public Scattering {
  double constant;
  AcousticScattering(const Material &m, const Band &b) : Scattering(m, b, 0) {
    double temperature = 1;
    constant =
        pow(material.acoustic_deformation_potential, 2) * consts::kB * temperature *
        (math::pi * pow(consts::hbar, 4) * material.density * pow(material.sound_velocity, 2));
  }
  virtual double rate(const Vec3 &p) const {
    return constant * band.acoustic_scattering_integral(band.energy(p));
  }
  virtual Vec3 scatter(const Vec3 &p, std::function<double()> random) const {
    return band.acoustic_scatter(band.energy(p), random);
  };
};

struct OpticalScattering : public Scattering {
  double constant;
  OpticalScattering(const Material &m, const Band &b, double energy)
      : Scattering(m, b, energy), constant(1) {}
  virtual double rate(const Vec3 &p) const {
    return constant * band.optical_scattering_integral(band.energy(p) - energy, p, energy);
  }
  virtual Vec3 scatter(const Vec3 &p, std::function<double()> random) const {
    return band.optical_scatter(band.energy(p) - energy, random);
  };
};

/**
  Dumper is used to collect data from simulation.
  Method `dump` is called at the end of each step of simulation.
*/
class Dumper {
public:
  virtual void dump(int particle, int step, const Particle &p) = 0;
  virtual ~Dumper(){};
};

void simulate(const std::vector<Scattering *> mechanisms,
              const Material &material,
              double temperature,
              std::function<Vec3(double, const Particle &)> force, // force(t, state) for rhs
              double time_step,
              double all_time,
              size_t ensemble_size,
              Dumper &dumper);

// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable
