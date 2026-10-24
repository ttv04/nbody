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

// N light bodies on circular orbits around a heavy central point mass.
// Body 0 is the central mass, the remaining N sit in an annulus.
class PointMassRing2D: public Generator {
public:
    PointMassRing2D(
        size_t N,
        double central_mass,
        double particle_mass,
        double inner_radius,
        double outer_radius,
        double G,
        double center_x = 0.0,
        double center_y = 0.0
    );
    State getRandomState();
    State getState(unsigned int seed);

private:
    size_t N_;
    double central_mass_;
    double particle_mass_;
    double inner_radius_;
    double outer_radius_;
    double G_;
    double center_x_;
    double center_y_;
};

// Dark-matter halo (truncated NFW), stellar disk (exponential), and
// central bulge (truncated Hernquist), in equilibrium in their shared potential.
struct GalaxyProfile {
    std::size_t halo_count = 0;
    double halo_mass = 0.0;
    double halo_scale_radius = 1.0;
    double halo_virial_radius = 10.0;

    std::size_t disk_count = 0;
    double disk_mass = 0.0;
    double disk_scale_radius = 1.0;
    double disk_scale_height = 0.1;
    double disk_outer_radius = 5.0;

    std::size_t bulge_count = 0;
    double bulge_mass = 0.0;
    double bulge_scale_radius = 1.0;
    double bulge_outer_radius = 10.0;

    Config halo_style{80, 90, 120, 255, 1.0f};
    Config disk_style{180, 200, 255, 255, 1.8f};
    Config bulge_style{255, 220, 170, 255, 2.2f};

    // Toomre Q of the stellar disk. Above 1 resists axisymmetric collapse.
    double toomre_q = 2.0;
    // Radians. Inclination tips the disk about x after the in-plane rotation.
    double inclination = 0.0;
    double position_angle = 0.0;
};

class CompositeGalaxy3D: public Generator {
public:
    CompositeGalaxy3D(
        GalaxyProfile profile,
        double gravitational_constant,
        double softening,
        double center_x = 0.0,
        double center_y = 0.0,
        double center_z = 0.0
    );
    State getRandomState();
    State getState(unsigned int seed);

private:
    GalaxyProfile profile_;
    double gravitational_constant_;
    double softening_;
    double center_x_;
    double center_y_;
    double center_z_;
};
