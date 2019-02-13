#pragma once

#include <ostream>
#include <vector>

#include <vec3.hh>

namespace probable {

struct Material {
  double mass;

  double energy(const Vec3 &p) const { return p.dot(p) / (2 * mass); }
  Vec3 velocity(const Vec3 &p) const { return p / mass; }
  Vec3 create_particle(double temperature) const;
};

struct Scattering {
  const Material &m;
  const double energy;
  Scattering(const Material &m, double e) : m(m), energy(e) {}
  virtual double rate(const Vec3 &p) const = 0;
  Vec3 scatter(const Vec3 &p) const;
  virtual ~Scattering() {}
};

struct Results {
  std::vector<Vec3> average_velocity;
  std::vector<double> average_power;
  std::vector<std::vector<size_t>> scatterings;
  Results() {}
  Results(size_t size)
      : average_velocity(size, Vec3()), average_power(size, 0), scatterings(size) {}

  friend std::ostream &operator<<(std::ostream &s, const Results &r);
};

Results simulate(const Material &material,
                 const std::vector<Scattering *> mechanisms,
                 double temperature,
                 const Vec3 &electric_field,
                 const Vec3 &magnetic_field,
                 double time_step,
                 double all_time,
                 size_t ansemble_size);

// Thread-safe Mersenne twister-based rng
double uniform();
} // namespace probable

namespace units {
// energy
const double eV = 1;
const double J = 1 / 1.6e-19;

// length
const double um = 1;
const double cm = um * 1e4;
const double m = um * 1e6;

// time
const double ps = 1;
const double s = ps * 1e12;

// mass
const double kg = J * s * s / m / m;
const double g = 1e-3 * kg;

// potential
const double V = 1;

// electric charge
const double C = J / V;

// temperature
const double K = 1.38e-23 / 1.6e-19;
} // namespace units

namespace consts {
const double hbar = 1.05e-34 * units::J * units::s;
const double e = 1.6e-19 * units::C;
const double c = 3e8 * units::m / units::s;
const double me = 9.1e-31 * units::kg;
const double kB = 1.38e-23 * units::J / units::K;
const double eps0 = 8.85e-12 * units::C / units::V / units::m;
} // namespace consts

namespace math {
const double e = exp(1);
const double pi = acos(-1);
} // namespace math
