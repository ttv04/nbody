#pragma once
#include <cstddef>
#include <initializer_list>
#include <vector>

class Config {
public:
    unsigned char r = 230;
    unsigned char g = 41;
    unsigned char b = 55;
    unsigned char a = 255;
    float radius = 5.0f;

    Config() = default;
    Config(
        unsigned char r,
        unsigned char g,
        unsigned char b,
        unsigned char a,
        float radius
    );
};

inline Config default_body_config{230, 41, 55, 255, 5.0f};

class Body {
public:
    std::vector<double> r;
    std::vector<double> v;
    double m;
    Config config;

    Body(
        std::initializer_list<double> r,
        std::initializer_list<double> v,
        double m,
        Config config = default_body_config
    );
    Body(
        std::size_t dimensions,
        double mass,
        Config config = default_body_config
    );

    std::size_t dimensions() const;
};

class State {
public:
    std::size_t dim;
    std::vector<Body> bodies;

    State(std::initializer_list<Body> bodies);
    State(std::size_t body_count, std::size_t dimensions, double mass = 1.0);
    State(const State& other);

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