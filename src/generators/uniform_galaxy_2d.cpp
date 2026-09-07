#include "nbody/generators.hpp"
#include <cstdlib>
#include <random>
#include <cmath>
#include <numbers>

/*
class UniformMassGalaxy2D: public Generator {
    public:
        UniformMassGalaxy2D(size_t N, double radius, double mass);
        State& getRandomState();
        State& getState(unsigned int seed);
    
    private:
        size_t N_;
        double radius_;
        double mass_;
};
*/

UniformMassGalaxy2D::UniformMassGalaxy2D(
    size_t N, 
    double radius,
    double mass,
    double center_x,
    double center_y,
    double G,
    double inner_radius
)
    : N_(N),
      radius_(radius),
      mass_(mass),
      center_x_(center_x),
      center_y_(center_y),
      G_(G),
      inner_radius_(inner_radius) {}

State UniformMassGalaxy2D::getRandomState() {
    std::random_device random;
    return getState(random());
}

State UniformMassGalaxy2D::getState(unsigned int seed) {
    std::mt19937 rng(seed);

    State state(N_, 2, mass_);

    double r_in2 = inner_radius_ * inner_radius_;
    double r_out2 = radius_ * radius_;
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::uniform_real_distribution<double> angle(0.0, 2 * std::numbers::pi);

    for (Body& body : state.bodies) {
        double r = sqrt(unit(rng) * (r_out2 - r_in2) + r_in2);
        double theta = angle(rng);

        body.r[0] = center_x_ + r * cos(theta);
        body.r[1] = center_y_ + r * sin(theta);

        double enclosed = mass_;
        if (r_out2 != r_in2) {
            enclosed = mass_ * (r * r - r_in2) / (r_out2 - r_in2);
        }

        double v_coeff = 0.0;
        if (r > 0.0) {
            v_coeff = sqrt(G_ * enclosed / r);
        }

        body.v[0] = -v_coeff * sin(theta);
        body.v[1] = v_coeff * cos(theta);

        body.m = mass_ / N_;
    }

    return state;
};
