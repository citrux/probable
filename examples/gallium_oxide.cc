#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <probable.hh>
#include <units.hh>

using namespace probable;
const double density = 5.88 * units::g / pow(units::cm, 3);     // density
const double sound_velocity = 6.8e3 * units::m / units::s;      // sound velocity
const double acoustic_deformation_potential = 16.6 * units::eV; // acoustic deformation potential
const double nonpolar_optical_deformation_potential =
    8.5e8 * units::eV / units::cm; // acoustic deformation potential
const double eps_inf = 4.21;       // high-frequency dielectric constant
const double eps_static = 11.4;    // static dielectric constant;

template <typename T> T parse(const std::string &s) {
  std::stringstream ss(s);
  T result;
  ss >> result;
  return result;
}

template <typename T> std::ostream &operator<<(std::ostream &s, std::vector<T> t) {
  s << "[";
  for (std::size_t i = 0; i < t.size(); i++) {
    s << t[i] << (i == t.size() - 1 ? "" : ", ");
  }
  return s << "]" << std::endl;
}

template <typename T> T sum(std::vector<T> t) {
  T s{};
  for (auto &v : t) {
    s += v;
  }
  return s;
}

class MyDumper: public Dumper {
  std::vector<Vec3> average_velocities;

public:
  MyDumper(size_t n) : average_velocities(n) {}

  void dump(int i, int step, const Particle &particle) {
    average_velocities[i] += (particle.band.velocity(particle.p) - average_velocities[i]) / (step + 1);
  }

  Vec3 mean() {
    Vec3 result;
    size_t n = average_velocities.size();
    for (size_t i = 0; i < n; ++i) {
      result += (average_velocities[i] - result) / (i + 1);
    }
    return result;
  }

  Vec3 std() {
    Vec3 result;
    size_t n = average_velocities.size();
    for (size_t i = 0; i < n; ++i) {
      result += (average_velocities[i] * average_velocities[i] - result) / (i + 1);
    }
    Vec3 m = mean();
    return (result - m * m).sqrt();
  }
};

// FIXME: global variables
// maybe make Simulation class with pure virtual Vec3 force() method
Vec3 electric_field;
Vec3 magnetic_field;
Vec3 force(double t, const Particle& p) {
  return p.charge * (electric_field + p.band.velocity(p.p).cross(magnetic_field));
}

int main(int argc, char const *argv[]) {
  if (argc != 4) {
    std::cout << "Invalid number of arguments\n";
    std::cout << "Usage: " << argv[0] << " <ensemble size> <temperature> <electric field>\n";
    return 1;
  }

  size_t ensemble_size = parse<int>(argv[1]);
  double temperature = parse<double>(argv[2]) * units::K;
  electric_field = {parse<double>(argv[3]) * units::V / units::m, 0, 0};
  magnetic_field = {0, 0, 0};

  Material gallium_oxide;
  gallium_oxide.bands = {new ParabolicBand{false, 0.29 * consts::me}};
  Band& conductive_band = *(gallium_oxide.bands[0]);
  std::vector<Scattering *> scattering_mechanisms{};
      // new AcousticScattering(gallium_oxide, conductive_band, temperature),
      // new ImpurityScattering(gallium_oxide, conductive_band, temperature),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 25e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 25e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 29e-3 * units::eV, 1.6e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 29e-3 * units::eV, 1.6e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 35e-3 * units::eV, 0.5e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 35e-3 * units::eV, 0.5e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 43e-3 * units::eV, 1.15e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 43e-3 * units::eV, 1.15e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 50e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 50e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 70e-3 * units::eV, 4.5e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 70e-3 * units::eV, 4.5e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 80e-3 * units::eV, 3.1e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 80e-3 * units::eV, 3.1e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 90e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 90e-3 * units::eV, 3.2e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 14e-3 * units::eV, 0.15e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 14e-3 * units::eV, 0.15e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 37e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 37e-3 * units::eV, 2.0e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 60e-3 * units::eV, 3.8e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 60e-3 * units::eV, 3.8e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalAbsorptionScattering(
      //     gallium_oxide, conductive_band, temperature, 81e-3 * units::eV, 3.0e2 * units::eV / pow(units::nm, 2)),
      // new PolarOpticalEmissionScattering(
      //     gallium_oxide, conductive_band, temperature, 81e-3 * units::eV, 3.0e2 * units::eV / pow(units::nm, 2))};

  double time_step = 1e-16 * units::s;
  double all_time = 1e-11 * units::s;

  // print info on stdout
  std::cout << "Ensemble size:    " << ensemble_size << "\n";
  std::cout << "Time step:        " << time_step / units::s << " s\n";
  std::cout << "Simulation time:  " << all_time / units::s << " s\n";
  std::cout << "Temperature:      " << temperature / units::K << " K\n";
  std::cout << "Electric field:   " << electric_field / units::V * units::m << " V/m\n";
  std::cout << "Magnetic field:   " << magnetic_field / units::T << " T\n";
  std::cout << "Scattering mechanisms:\n";
  // for (size_t i = 0; i < scattering_mechanisms.size(); ++i) {
  //   std::cout << i + 1 << ": " << *scattering_mechanisms[i] << '\n';
  // }

  MyDumper dumper(ensemble_size);
  simulate(scattering_mechanisms,
           gallium_oxide,
                          temperature,
                          force,
                          time_step,
                          all_time,
                          ensemble_size,
                          dumper);
  Vec3 average_velocity;
  Vec3 average_velocity2;
  //  std::vector<double> scattering_rates(scattering_mechanisms.size(), 0);
  std::cout << "Average velocity: " << dumper.mean() << '\n';
  std::cout << "             std: " << dumper.std() << '\n';
  // std::cout << "Scattering rates: " << scattering_rates << '\n';
  // std::cout << "           Total: " << sum(scattering_rates) << '\n';
  return 0;
}
