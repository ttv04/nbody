#include "nbody/state.hpp"

class Generator {
    virtual State getRandomState() = 0;
    virtual State getState(unsigned int seed) = 0;
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