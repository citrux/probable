#include <omp.h>
#include <random>

#include <probable.hh>

namespace probable {

Vec3 Material::create_particle(double temperature) const {
  double p_max = 5 * sqrt(2 * mass * consts::kB * temperature);
  while (true) {
    double prob = uniform();
    double p1 = p_max * cbrt(uniform());
    if (prob < exp(-p1 * p1 / (2 * mass * consts::kB * temperature))) {
      double cos_theta = 1 - 2 * uniform();
      double sin_theta = sqrt(1 - cos_theta * cos_theta);
      double phi = 2 * math::pi * uniform();
      Vec3 p = {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
      p *= p1;
      return p;
    }
  }
}

Vec3 Scattering::scatter(const Vec3 &p) const {
  double e = m.energy(p) - energy;
  double r = sqrt(2 * m.mass * e);
  double cos_theta = 1 - 2 * uniform();
  double sin_theta = sqrt(1 - cos_theta * cos_theta);
  double phi = 2 * math::pi * uniform();
  return {r * sin_theta * cos(phi), r * sin_theta * sin(phi), r * cos_theta};
}

std::ostream &operator<<(std::ostream &s, const Results &r) {
  size_t size = r.average_velocity.size();
  for (size_t i = 0; i < size; ++i) {
    s << r.average_velocity[i].x / units::m * units::s << " "
      << r.average_velocity[i].y / units::m * units::s << " "
      << r.average_velocity[i].z / units::m * units::s << " "
      << r.average_power[i] / units::J * units::s;
    for (size_t v : r.scatterings[i]) {
      s << " " << v;
    }
    s << "\n";
  }
  return s;
}

Results simulate(const Material &material,
                 const std::vector<Scattering *> mechanisms, double temperature,
                 const Vec3 &electric_field, const Vec3 &magnetic_field,
                 double time_step, double all_time, size_t ansemble_size) {
  Results results(ansemble_size);
  size_t steps = all_time / time_step + 1;
#pragma omp parallel for
  for (size_t i = 0; i < ansemble_size; ++i) {
    Vec3 &average_velocity = results.average_velocity[i];
    double &average_power = results.average_power[i];
    std::vector<size_t> &scatterings = results.scatterings[i];

    average_velocity = {0, 0, 0};
    average_power = 0;
    scatterings.assign(mechanisms.size(), 0);
    Vec3 p = material.create_particle(temperature);
    std::vector<double> free_flight(mechanisms.size(), 0);
    for (double &l : free_flight) {
      l = -log(uniform());
    }
    for (size_t j = 0; j < steps; ++j) {
      Vec3 v = material.velocity(p);
      average_velocity += (v - average_velocity) / (j + 1);
      average_power +=
          (v.dot(electric_field) - results.average_power[i]) / (j + 1);
      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p) * time_step;
        if (free_flight[k] < 0) {
          p = mechanisms[k]->scatter(p);
          scatterings[k] += 1;
          free_flight[k] = -log(uniform());
          break;
        }
      }
      p += -consts::e * (electric_field + v.cross(magnetic_field)) * time_step;
    }
  }
  return results;
}

double uniform() {
  thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
}
} // namespace probable
