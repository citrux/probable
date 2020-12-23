#pragma once

#include <cmath>

namespace units {
// energy
const double eV = 1;
const double J = 1 / 1.6e-19;

// length
const double um = 1;
const double cm = um * 1e4;
const double m = um * 1e6;
const double nm = um * 1e-3;

// time
const double ps = 1;
const double s = ps * 1e12;

// mass
const double kg = J * s * s / m / m;
const double g = 1e-3 * kg;

// potential
const double V = 1;

// magnetic field
const double T = 3e8 * (V / m) / (m / s);

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

