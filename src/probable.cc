#include <iostream>
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

void Results::append(uint32_t n, double t, const Vec3 &p, const Vec3 &v, double e, size_t s) {
  if (not(flags & DumpFlags::on_scatterings) or s) {
    if (flags & DumpFlags::scattering) {
      scatterings.push_back(s);
    }
    if (flags & DumpFlags::number) {
      ns.push_back(n);
    }
    if (flags & DumpFlags::time) {
      ts.push_back(t);
    }
    if (flags & DumpFlags::momentum) {
      momentums.push_back(p);
    }
    if (flags & DumpFlags::velocity) {
      velocities.push_back(v);
    }
    if (flags & DumpFlags::energy) {
      energies.push_back(e);
    }
    size += 1;
  }
}

std::ostream &operator<<(std::ostream &s, const Results &r) {
  for (size_t i = 0; i < r.size; ++i) {
    if (r.flags & DumpFlags::number) {
      s << r.ns[i] << " ";
    }
    if (r.flags & DumpFlags::time) {
      s << r.ts[i] / units::s << " ";
    }
    if (r.flags & DumpFlags::momentum) {
      s << r.momentums[i].x << " ";
      s << r.momentums[i].y << " ";
      s << r.momentums[i].z << " ";
    }
    if (r.flags & DumpFlags::velocity) {
      s << r.velocities[i].x / units::m * units::s << " ";
      s << r.velocities[i].y / units::m * units::s << " ";
      s << r.velocities[i].z / units::m * units::s << " ";
    }
    if (r.flags & DumpFlags::energy) {
      s << r.energies[i] / units::eV << " ";
    }
    if (r.flags & DumpFlags::scattering) {
      s << r.scatterings[i];
    }
    s << "\n";
  }
  return s;
}

std::vector<Results> simulate(const Material &material,
                              const std::vector<Scattering *> mechanisms,
                              double temperature,
                              const Vec3 &electric_field,
                              const Vec3 &magnetic_field,
                              double time_step,
                              double all_time,
                              size_t ansemble_size,
                              DumpFlags flags) {
  std::vector<Results> results(ansemble_size);
  size_t steps = all_time / time_step + 1;
  size_t alloc = (flags & DumpFlags::on_scatterings) ? steps / 10 : steps;
#pragma omp parallel for
  for (size_t i = 0; i < ansemble_size; ++i) {
    Results &result = results[i];
    result = Results(alloc, flags);
    Vec3 p = material.create_particle(temperature);
    std::vector<double> free_flight(mechanisms.size(), 0);
    for (double &l : free_flight) {
      l = -log(uniform());
    }

    for (size_t j = 0; j < steps; ++j) {
      Vec3 p_ = p;
      Vec3 v = material.velocity(p_);
      double e = material.energy(p_);
      size_t scattering_mechanism = 0; // means no scattering
      for (size_t k = 0; k < mechanisms.size(); ++k) {
        free_flight[k] -= mechanisms[k]->rate(p_) * time_step;
        if (free_flight[k] < 0) {
          p = mechanisms[k]->scatter(p_);
          free_flight[k] = -log(uniform());
          scattering_mechanism = k + 1; // enumerate mechanisms from 1
          break;
        }
      }
      result.append(j, j * time_step, p_, v, e, scattering_mechanism);
      if (not scattering_mechanism) {
        p += -consts::e * (electric_field + v.cross(magnetic_field)) * time_step;
      }
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
