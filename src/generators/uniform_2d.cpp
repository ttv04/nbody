#include "nbody/generators.hpp"
#include <random>

UniformMass2D::UniformMass2D(
    size_t N,
    double starting_box_width,
    double mass_per_particle,
    double center_x,
    double center_y
)
    : N_(N),
      starting_box_width_(starting_box_width),
      mass_per_particle_(mass_per_particle),
      center_x_(center_x),
      center_y_(center_y) {}

State UniformMass2D::getRandomState() {
    std::random_device random;
    return getState(random());
}

State UniformMass2D::getState(unsigned int seed) {
    std::mt19937 rng(seed);

    State state(N_, 2, mass_per_particle_);
    const double half = 0.5 * starting_box_width_;
    std::uniform_real_distribution<double> x(center_x_ - half, center_x_ + half);
    std::uniform_real_distribution<double> y(center_y_ - half, center_y_ + half);

    for (Body& body : state.bodies) {
        body.r[0] = x(rng);
        body.r[1] = y(rng);
        body.v[0] = 0.0;
        body.v[1] = 0.0;
        body.m = mass_per_particle_;
    }

    return state;
}
