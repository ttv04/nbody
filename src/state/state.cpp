#include "nbody/state.hpp"

#include <algorithm>

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