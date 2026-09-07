#pragma once
#include <cstddef>
#include <initializer_list>
#include <vector>

class Body {
public:
    std::vector<double> r;
    std::vector<double> v;
    double m;

    Body(std::initializer_list<double> r, std::initializer_list<double> v, double m);
    Body(std::size_t dimensions, double mass);

    std::size_t dimensions() const;
};

class State {
public:
    std::size_t dim;
    std::vector<Body> bodies;

    State(std::initializer_list<Body> bodies);
    State(std::size_t body_count, std::size_t dimensions, double mass = 1.0);

    std::size_t size() const;
};

class AccelerationField {
public:
    AccelerationField() = default;
    AccelerationField(std::size_t body_count, std::size_t dimensions);

    void resize(std::size_t body_count, std::size_t dimensions);
    void reset();

    std::size_t body_count() const;
    std::size_t dimensions() const;

    double& operator()(std::size_t body, std::size_t dimension);
    const double& operator()(std::size_t body, std::size_t dimension) const;

private:
    std::size_t body_count_{0};
    std::size_t dimensions_{0};
    std::vector<double> values_;
};