#include "nbody/state.hpp"

#include <algorithm>
#include <cmath>

Config::Config(
    unsigned char r,
    unsigned char g,
    unsigned char b,
    unsigned char a,
    float radius
)
    : r(r), g(g), b(b), a(a), radius(radius) {}

Body::Body(
    std::initializer_list<double> position,
    std::initializer_list<double> velocity,
    double mass,
    Config config
)
    : r(position), v(velocity), m(mass), config(config) {}

Body::Body(std::size_t dimensions, double mass, Config config)
    : r(dimensions, 0.0), v(dimensions, 0.0), m(mass), config(config) {}

std::size_t Body::dimensions() const {
    return r.size();
}

State::State(std::initializer_list<Body> initial_bodies)
    : dim(initial_bodies.size() == 0 ? 0 : initial_bodies.begin()->dimensions()),
      bodies(initial_bodies) {}

State::State(std::size_t body_count, std::size_t dimensions, double mass)
    : dim(dimensions), bodies(body_count, Body(dimensions, mass)) {}

State::State(const State& other)
    : dim(other.dim),
      bodies(other.bodies) {}

double State::total_potential_energy(
    double gravitational_constant,
    double softening
) const {
    double pe = 0.0;
    for (std::size_t i = 0; i < bodies.size(); i++) {
        for (std::size_t j = i + 1; j < bodies.size(); j++) {
            const Body& a = bodies[i];
            const Body& b = bodies[j];
            double r2 = softening * softening;
            for (std::size_t d = 0; d < dim; d++) {
                const double displacement = b.r[d] - a.r[d];
                r2 += displacement * displacement;
            }
            if (r2 == 0.0) {
                continue;
            }
            pe -= gravitational_constant * a.m * b.m / std::sqrt(r2);
        }
    }
    return pe;
}

double State::total_kinetic_energy() const {
    double ke = 0.0;
    for (std::size_t i = 0; i < bodies.size(); i++) {
        const Body& body = bodies[i];
        double v2 = 0.0;
        for (std::size_t d = 0; d < dim; d++) {
            v2 += body.v[d] * body.v[d];
        }
        ke += body.m * v2;
    }

    return 0.5 * ke;
}

double State::total_energy(
    double gravitational_constant,
    double softening
) const {
    return total_kinetic_energy()
        + total_potential_energy(gravitational_constant, softening);
}

std::size_t State::size() const {
    return bodies.size();
}

AccelerationField::AccelerationField(
    std::size_t body_count,
    std::size_t dimensions
) {
    resize(body_count, dimensions);
}

void AccelerationField::resize(
    std::size_t body_count,
    std::size_t dimensions
) {
    body_count_ = body_count;
    dimensions_ = dimensions;
    values_.assign(body_count * dimensions, 0.0);
}

void AccelerationField::reset() {
    std::fill(values_.begin(), values_.end(), 0.0);
}

std::size_t AccelerationField::body_count() const {
    return body_count_;
}

std::size_t AccelerationField::dimensions() const {
    return dimensions_;
}

double& AccelerationField::operator()(
    std::size_t body,
    std::size_t dimension
) {
    return values_[body * dimensions_ + dimension];
}

const double& AccelerationField::operator()(
    std::size_t body,
    std::size_t dimension
) const {
    return values_[body * dimensions_ + dimension];
}