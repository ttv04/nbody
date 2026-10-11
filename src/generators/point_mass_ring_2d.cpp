#include "nbody/generators.hpp"

#include <cmath>
#include <numbers>
#include <random>

PointMassRing2D::PointMassRing2D(
    size_t N,
    double central_mass,
    double particle_mass,
    double inner_radius,
    double outer_radius,
    double G,
    double center_x,
    double center_y
)
    : N_(N),
      central_mass_(central_mass),
      particle_mass_(particle_mass),
      inner_radius_(inner_radius),
      outer_radius_(outer_radius),
      G_(G),
      center_x_(center_x),
      center_y_(center_y) {}

State PointMassRing2D::getRandomState() {
    std::random_device random;
    return getState(random());
}

State PointMassRing2D::getState(unsigned int seed) {
    std::mt19937 rng(seed);

    State state(N_ + 1, 2, particle_mass_);

    Body& central = state.bodies[0];
    central.r[0] = center_x_;
    central.r[1] = center_y_;
    central.v[0] = 0.0;
    central.v[1] = 0.0;
    central.m = central_mass_;
    central.config = Config{255, 196, 92, 255, 5.0f};

    const double r_in2 = inner_radius_ * inner_radius_;
    const double r_out2 = outer_radius_ * outer_radius_;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> angle(0.0, 2 * std::numbers::pi);

    for (size_t i = 1; i < state.size(); i++) {
        // Uniform in area so the ring does not pile up at the inner edge.
        const double r = std::sqrt(unit(rng) * (r_out2 - r_in2) + r_in2);
        const double theta = angle(rng);

        // Circular orbit speed for the central mass alone; the ring is light.
        const double v = std::sqrt(G_ * central_mass_ / r);

        Body& body = state.bodies[i];
        body.r[0] = center_x_ + r * std::cos(theta);
        body.r[1] = center_y_ + r * std::sin(theta);
        body.v[0] = -v * std::sin(theta);
        body.v[1] = v * std::cos(theta);
        body.m = particle_mass_;
        body.config = Config{170, 205, 255, 255, 0.6f};
    }

    return state;
}
