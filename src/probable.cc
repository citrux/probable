#include <iostream>
#include <omp.h>
#include <random>

#include <probable.hh>
#include <units.hh>

namespace probable {

Particle Material::create_particle(double temperature) const {
  while (true) {
    for (auto band : bands) {
      Vec3 p;
      if (band->try_boltzmann_sample_momentum(temperature, p)) {
        return Particle{p, Vec3(), consts::e, *band};
      }
    }
  }
}

bool ParabolicBand::try_boltzmann_sample_momentum(double temperature, Vec3 &momentum) const {
  double p_max = 5 * sqrt(2 * mass * consts::kB * temperature);
  double prob = uniform();
  double p1 = p_max * cbrt(uniform());
  double cos_theta = 1 - 2 * uniform();
  double sin_theta = sqrt(1 - cos_theta * cos_theta);
  double phi = 2 * math::pi * uniform();
  Vec3 p = {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
  p *= p1;
  if (prob < exp(-energy(p) / (consts::kB * temperature))) {
    momentum = p;
    return true;
  }
  return false;
}

void simulate(const std::vector<Scattering *> mechanisms,
              const Material &material,
              double temperature,
              std::function<Vec3(double, const Particle &)> force, // force(t, state) for rhs
              double time_step,
              double all_time,
              size_t ensemble_size,
              Dumper &dumper) {
  size_t steps = all_time / time_step + 1;
#pragma omp parallel for
  for (size_t i = 0; i < ensemble_size; ++i) {
    Particle p = material.create_particle(temperature);
    std::vector<double> free_flight(mechanisms.size(), 0);
    for (double &l : free_flight) {
      l = -log(uniform());
    }

    for (size_t j = 0; j < steps; ++j) {
      Vec3 p_ = p.p;
      Vec3 v = p.band.velocity(p_);
      double e = p.band.energy(p_);
      size_t scattering_mechanism = 0; // means no scattering
      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p_) * time_step;
        if (free_flight[k] < 0) {
          p.p = mechanisms[k]->scatter(p_, uniform);
          free_flight[k] = -log(uniform());
          scattering_mechanism = k + 1; // enumerate mechanisms from 1
          break;
        }
      }
      dumper.dump(i, j, p);
      if (not scattering_mechanism) {
        p.p += force(j * time_step, p) * time_step;
      }
    }
  }
}

double uniform() {
  thread_local std::mt19937 generator(std::random_device{}());
  std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
}
} // namespace probable
