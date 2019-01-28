#include <iostream>
#include <sstream>

#include <probable.hh>

using namespace probable;
const double density = 5.88 * units::g / pow(units::cm, 3); // density           
const double sound_velocity = 6.8e3 * units::m / units::s; // sound velocity
const double acoustic_deformation_potential = 16.6 * units::eV; // acoustic deformation potential
const double nonpolar_optical_deformation_potential = 8.5e8 * units::eV / units::cm; // acoustic deformation potential
const double eps_inf = 4.21;    // high-frequency dielectric constant
const double eps_static = 11.4; // static dielectric constant;

struct AcousticScattering : public Scattering {
  AcousticScattering(const Material &m, double temperature) : Scattering(m, 0) {
    constant = pow(acoustic_deformation_potential, 2) * consts::kB *
               temperature * m.mass /
               (math::pi * pow(consts::hbar, 4) * density *
                pow(sound_velocity, 2));
  }
  double rate(const Vec3 &p) const { return constant * p.length(); }

private:
  double constant;
};

struct PolarOpticalAbsorptionScattering : public Scattering {
  PolarOpticalAbsorptionScattering(const Material &m, double temperature,
                                   double energy)
      : Scattering(m, -energy) {
    constant = pow(consts::e, 2) / (4 * math::pi * consts::eps0) * energy /
               (2 * math::pi * pow(consts::hbar, 2)) *
               (1 / eps_inf - 1 / eps_static) * 1 /
               (exp(energy / consts::kB / temperature) - 1);
  }
  double rate(const Vec3 &p) const {
    double v = m.velocity(p).length();
    if (v < 1e-20) {
      return constant * sqrt(m.mass / std::abs(2 * energy));
    }
    return constant / v * asinh(sqrt(m.energy(p) / std::abs(energy)));
  }

private:
  double constant;
};

struct PolarOpticalEmissionScattering : public Scattering {
  PolarOpticalEmissionScattering(const Material &m, double temperature,
                                 double energy)
      : Scattering(m, energy) {
    constant = pow(consts::e, 2) / (4 * math::pi * consts::eps0) * energy /
               (2 * math::pi * pow(consts::hbar, 2)) *
               (1 / eps_inf - 1 / eps_static) *
               (1 + 1 / (exp(energy / consts::kB / temperature) - 1));
  }
  double rate(const Vec3 &p) const {
    double e = m.energy(p);
    if (e < energy) {
      return 0;
    }
    double v = m.velocity(p).length();
    return constant / v * acosh(sqrt(e / std::abs(energy)));
  }

private:
  double constant;
};

template <typename T> T parse(const std::string &s) {
  std::stringstream ss(s);
  T result;
  ss >> result;
  return result;
}

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0]
              << " <ansemble size> <temperature> <electric field>\n";
    return 1;
  }

  size_t ansemble_size = parse<int>(argv[1]);
  double temperature = parse<double>(argv[2]) * units::K;
  Vec3 electric_field{parse<double>(argv[3]) * units::V / units::m, 0, 0};
  Vec3 magnetic_field{0, 0, 0};

  Material gallium_oxide{0.29 * consts::me}; // electron effective mass

  std::vector<Scattering *> scattering_mechanisms{
      new AcousticScattering(gallium_oxide, temperature),
      new PolarOpticalAbsorptionScattering(gallium_oxide, temperature,
                                           21e-3 * units::eV),
      new PolarOpticalEmissionScattering(gallium_oxide, temperature,
                                         21e-3 * units::eV),
      new PolarOpticalAbsorptionScattering(gallium_oxide, temperature,
                                           44e-3 * units::eV),
      new PolarOpticalEmissionScattering(gallium_oxide, temperature,
                                         44e-3 * units::eV)};

  double time_step = 1e-15 * units::s;
  double all_time = 1e-9 * units::s;
  auto results = simulate(gallium_oxide, scattering_mechanisms, temperature,
                          electric_field, magnetic_field, time_step, all_time,
                          ansemble_size);
  std::cout << results;
  return 0;
}