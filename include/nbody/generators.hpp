#pragma once

#include "nbody/state.hpp"

class Generator {
    virtual State getRandomState() = 0;
    virtual State getState(unsigned int seed) = 0;
};

class UniformMass2D: public Generator {
public:
    UniformMass2D(
        size_t N,
        double starting_box_width,
        double mass_per_particle,
        double center_x = 0.0,
        double center_y = 0.0
    );
    State getRandomState();
    State getState(unsigned int seed);
private:
    size_t N_;
    double starting_box_width_;
    double mass_per_particle_;
    double center_x_;
    double center_y_;
};

class UniformMass3D: public Generator {
public:
    UniformMass3D(
        size_t N,
        double starting_box_width,
        double mass_per_particle,
        double center_x = 0.0,
        double center_y = 0.0,
        double center_z = 0.0
    );
    State getRandomState();
    State getState(unsigned int seed);
private:
    size_t N_;
    double starting_box_width_;
    double mass_per_particle_;
    double center_x_;
    double center_y_;
    double center_z_;
};

class UniformMassGalaxy2D: public Generator {
public:
    UniformMassGalaxy2D(size_t N, double radius, double mass, double center_x, double center_y, double G, double inner_radius = 0.0);
    State getRandomState();
    State getState(unsigned int seed);

private:
    size_t N_;
    double radius_;
    double mass_;
    double center_x_;
    double center_y_;
    double G_;
    double inner_radius_;
};

class UniformMassGalaxy3D: public Generator {
public:
    UniformMassGalaxy3D(
        size_t N,
        double radius,
        double mass,
        double center_x,
        double center_y,
        double center_z,
        double G,
        double inner_radius = 0.0
    );
    State getRandomState();
    State getState(unsigned int seed);

private:
    size_t N_;
    double radius_;
    double mass_;
    double center_x_;
    double center_y_;
    double center_z_;
    double G_;
    double inner_radius_;
};
