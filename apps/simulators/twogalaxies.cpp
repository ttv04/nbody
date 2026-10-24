#include "nbody/core.hpp"
#include "nbody/data.hpp"
#include "nbody/generators.hpp"
#include "argparse/argparse.hpp"
#include "sim_options.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <vector>

namespace {

const std::size_t halo_count = 3000;
const std::size_t disk_count = 500;
const std::size_t bulge_count = 150;
const double halo_mass = 800.0;
const double disk_mass = 80.0;
const double bulge_mass = 20.0;
const double halo_scale_radius = 6.0;
const double halo_virial_radius = 60.0;
const double disk_scale_radius = 3.5;
const double disk_scale_height = 0.4;
const double disk_outer_radius = 17.5;
const double bulge_scale_radius = 0.8;
const double bulge_outer_radius = 12.0;
const double center_separation = 110.0;
const double impact_parameter = 5.0;

GalaxyProfile make_profile(
    Config halo_style,
    Config disk_style,
    Config bulge_style,
    double inclination,
    double position_angle
) {
    GalaxyProfile profile;
    profile.halo_count = halo_count;
    profile.halo_mass = halo_mass;
    profile.halo_scale_radius = halo_scale_radius;
    profile.halo_virial_radius = halo_virial_radius;
    profile.disk_count = disk_count;
    profile.disk_mass = disk_mass;
    profile.disk_scale_radius = disk_scale_radius;
    profile.disk_scale_height = disk_scale_height;
    profile.disk_outer_radius = disk_outer_radius;
    profile.bulge_count = bulge_count;
    profile.bulge_mass = bulge_mass;
    profile.bulge_scale_radius = bulge_scale_radius;
    profile.bulge_outer_radius = bulge_outer_radius;
    profile.halo_style = halo_style;
    profile.disk_style = disk_style;
    profile.bulge_style = bulge_style;
    profile.toomre_q = 2.0;
    profile.inclination = inclination;
    profile.position_angle = position_angle;
    return profile;
}

State make_galaxy(
    const GalaxyProfile& profile,
    double center_x,
    double center_y,
    double center_z,
    double bulk_vx,
    double bulk_vy,
    double bulk_vz,
    double gravitational_constant,
    double softening,
    unsigned seed
) {
    CompositeGalaxy3D galaxy(
        profile,
        gravitational_constant,
        softening,
        center_x,
        center_y,
        center_z
    );
    State state = galaxy.getState(seed);
    for (Body& body : state.bodies) {
        body.v[0] += bulk_vx;
        body.v[1] += bulk_vy;
        body.v[2] += bulk_vz;
    }
    return state;
}

State combine(State first, State second) {
    std::vector<Body> bodies = std::move(first.bodies);
    bodies.insert(
        bodies.end(),
        std::make_move_iterator(second.bodies.begin()),
        std::make_move_iterator(second.bodies.end())
    );
    return State(std::move(bodies));
}

}  // namespace

int main(int argc, const char* argv[]) {
    argparse::ArgumentParser program("two_galaxy");

    add_sim_args(program);

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    SimOptions o = sim_options_from(program);

    // Softening is near the disk thickness. The step is short compared
    // with an orbit in the stellar disk.
    if (!program.is_used("--dt")) {
        o.dt = 0.1;
    }
    if (!program.is_used("--soft")) {
        o.soft = 0.35;
    }
    if (!program.is_used("--solver")) {
        o.solver_name = "barnes-hut";
    }
    if (o.solver_name == "barnes-hut") {
        o.solver = std::make_unique<BarnesHutSolver>(
            o.soft,
            o.g,
            o.theta,
            o.max_leaf
        );
    }
    else {
        o.solver = std::make_unique<NaiveSolver>(o.soft, o.g);
    }

    const double galaxy_mass_total = halo_mass + disk_mass + bulge_mass;
    const double escape_speed = std::sqrt(
        2.0 * o.g * 2.0 * galaxy_mass_total / center_separation
    );
    const double approach_speed = 0.10 * escape_speed;

    const GalaxyProfile first_profile = make_profile(
        Config{70, 95, 150, 255, 1.0f},
        Config{110, 180, 255, 255, 1.8f},
        Config{255, 230, 190, 255, 2.2f},
        0.45,
        0.2
    );
    const GalaxyProfile second_profile = make_profile(
        Config{140, 85, 70, 255, 1.0f},
        Config{255, 160, 70, 255, 1.8f},
        Config{255, 214, 170, 255, 2.2f},
        1.05,
        0.8
    );

    State first = make_galaxy(
        first_profile,
        -0.5 * center_separation,
        0.5 * impact_parameter,
        0.0,
        approach_speed,
        0.0,
        0.0,
        o.g,
        o.soft,
        o.seed
    );
    State second = make_galaxy(
        second_profile,
        0.5 * center_separation,
        -0.5 * impact_parameter,
        0.0,
        -approach_speed,
        0.0,
        0.0,
        o.g,
        o.soft,
        o.seed + 1
    );
    const State init_state = combine(std::move(first), std::move(second));

    std::vector<double> masses;
    std::vector<Config> styles;
    masses.reserve(init_state.size());
    styles.reserve(init_state.size());
    for (const Body& body : init_state.bodies) {
        masses.push_back(body.m);
        styles.push_back(body.config);
    }

    Simulation sim(*o.integrator, *o.solver, init_state);

    int steps = (int)(o.t_start / o.dt);
    for (int i = 0; i < steps; i++) {
        sim.step(o.dt);
    }

    Result result{
        "Two Galaxies",
        o.integrator_name,
        o.solver_name,
        o.theta,
        o.max_leaf,
        init_state.size(),
        init_state.dim,
        o.g,
        o.soft,
        o.t_start,
        o.t_end,
        o.dt,
        static_cast<double>(o.seed),
        masses,
        styles,
        std::vector<double>(),
        std::vector<double>(),
        std::vector<std::vector<std::vector<double>>>(),
        std::vector<std::vector<std::vector<double>>>(),
    };

    steps = (int)((o.t_end - o.t_start) / o.dt);
    for (int i = 0; i < steps; i++) {
        sim.step(o.dt);
        const State& state = sim.state();

        std::vector<std::vector<double>> displacement;
        std::vector<std::vector<double>> velocity;
        displacement.reserve(state.size());
        velocity.reserve(state.size());

        for (std::size_t body_index = 0; body_index < state.size(); body_index++) {
            const Body& body = state.bodies[body_index];
            displacement.push_back(body.r);
            velocity.push_back(body.v);
        }

        result.t.push_back(o.t_start + o.dt * (i + 1));
        result.energy.push_back(state.total_energy(o.g, o.soft));
        result.displacement.push_back(displacement);
        result.velocity.push_back(velocity);
    }

    save_result(result, o.output);
}
