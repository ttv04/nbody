#include "nbody/generators.hpp"
#include <cmath>
#include <numbers>
#include <random>

UniformMassGalaxy3D::UniformMassGalaxy3D(
    size_t N,
    double radius,
    double mass,
    double center_x,
    double center_y,
    double center_z,
    double G,
    double inner_radius
)
    : N_(N),
      radius_(radius),
      mass_(mass),
      center_x_(center_x),
      center_y_(center_y),
      center_z_(center_z),
      G_(G),
      inner_radius_(inner_radius) {}

State UniformMassGalaxy3D::getRandomState() {
    std::random_device random;
    return getState(random());
}

State UniformMassGalaxy3D::getState(unsigned int seed) {
    std::mt19937 rng(seed);

    State state(N_, 3, mass_);

    const double r_in3 = inner_radius_ * inner_radius_ * inner_radius_;
    const double r_out3 = radius_ * radius_ * radius_;
    std::uniform_real_distribution<double> unit(0.0, 1.0);

    for (Body& body : state.bodies) {
        const double r = std::cbrt(unit(rng) * (r_out3 - r_in3) + r_in3);
        const double cos_theta = 2.0 * unit(rng) - 1.0;
        const double sin_theta = std::sqrt(std::max(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = 2.0 * std::numbers::pi * unit(rng);

        const double dx = r * sin_theta * std::cos(phi);
        const double dy = r * sin_theta * std::sin(phi);
        const double dz = r * cos_theta;

        body.r[0] = center_x_ + dx;
        body.r[1] = center_y_ + dy;
        body.r[2] = center_z_ + dz;

        double enclosed = mass_;
        if (r_out3 != r_in3) {
            enclosed = mass_ * (r * r * r - r_in3) / (r_out3 - r_in3);
        }

        double v_coeff = 0.0;
        if (r > 0.0) {
            v_coeff = std::sqrt(G_ * enclosed / r);
        }

        double tx = -dy;
        double ty = dx;
        double tz = 0.0;
        double tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
        if (tlen < 1e-12) {
            tx = 0.0;
            ty = -dz;
            tz = dy;
            tlen = std::sqrt(tx * tx + ty * ty + tz * tz);
        }
        if (tlen > 0.0) {
            tx /= tlen;
            ty /= tlen;
            tz /= tlen;
        }

        body.v[0] = v_coeff * tx;
        body.v[1] = v_coeff * ty;
        body.v[2] = v_coeff * tz;
        body.m = mass_ / N_;
    }

    return state;
}
