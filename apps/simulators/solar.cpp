#include "nbody/core.hpp"
#include "nbody/data.hpp"
#include "nbody/system.hpp"
#include "argparse/argparse.hpp"
#include "sim_options.hpp"

#include <exception>
#include <iostream>
#include <memory>

// CODATA 2022 Newtonian constant of gravitation, m^3 kg^-1 s^-2.
const double gravitational_constant = 6.67430e-11;

// NASA planetary fact sheet masses for the Sun, Earth, and Moon, kg.
const double sun_mass = 1.9885e30;
const double earth_mass = 5.9722e24;
const double moon_mass = 7.342e22;

// JPL planetary physical parameters, kg.
const double mercury_mass = 0.330103e24;
const double venus_mass = 4.86731e24;
const double mars_mass = 0.641691e24;
const double jupiter_mass = 1898.125e24;
const double saturn_mass = 568.317e24;
const double uranus_mass = 86.8099e24;
const double neptune_mass = 102.4092e24;

// Exact SI astronomical unit. The Earth starts this far from the Sun.
const double earth_sun_distance = 149597870700.0;

// JPL Keplerian semimajor axes at J2000, in AU, converted with the exact AU.
const double mercury_distance = 0.38709927 * earth_sun_distance;
const double venus_distance = 0.72333566 * earth_sun_distance;
const double mars_distance = 1.52371034 * earth_sun_distance;
const double jupiter_distance = 5.20288700 * earth_sun_distance;
const double saturn_distance = 9.53667594 * earth_sun_distance;
const double uranus_distance = 19.18916464 * earth_sun_distance;
const double neptune_distance = 30.06992276 * earth_sun_distance;

// NASA mean Earth-Moon distance, 384400 km.
const double earth_moon_distance = 3.844e8;

// Mean radii, used only as the drawn size. m.
const Config sun_style{255, 196, 92, 255, 6.957e8f};
const Config mercury_style{176, 168, 160, 255, 2.4394e6f};
const Config venus_style{214, 178, 92, 255, 6.0518e6f};
const Config earth_style{70, 130, 220, 255, 6.371e6f};
const Config moon_style{186, 186, 196, 255, 1.737e6f};
const Config mars_style{196, 92, 58, 255, 3.3895e6f};
const Config jupiter_style{206, 158, 112, 255, 6.9911e7f};
const Config saturn_style{224, 196, 140, 255, 5.8232e7f};
const Config uranus_style{148, 214, 220, 255, 2.5362e7f};
const Config neptune_style{62, 96, 214, 255, 2.4622e7f};

namespace {

// JPL publishes GM more precisely than mass. Mass is GM / G.
double mass_from_gm(double gm_km3_s2) {
    return gm_km3_s2 * 1.0e9 / gravitational_constant;
}

void add_moon(
    System& planet,
    double gm_km3_s2,
    double semimajor_axis_km,
    double theta,
    Config style,
    bool retrograde = false
) {
    planet.add(
        Body(2, mass_from_gm(gm_km3_s2), style),
        semimajor_axis_km * 1000.0,
        theta,
        retrograde
    );
}

}  // namespace

int main(int argc, const char* argv[]) {
    argparse::ArgumentParser program("solar");

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

    // Shared defaults are for the toy demos (G = 1, dt = 0.1 s, soft = 0.5).
    // This run is in seconds and meters, so keep an hour-long step and no
    // softening unless the command line set them. G is the CODATA value.
    if (!program.is_used("--dt")) {
        o.dt = 3600.0;
    }
    if (!program.is_used("--soft")) {
        o.soft = 0.0;
    }
    o.g = gravitational_constant;
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

    System earth_system(Body(2, earth_mass, earth_style), o.g);
    earth_system.add(Body(2, moon_mass, moon_style), earth_moon_distance, 0.0);

    // JPL planetary satellite GM (km^3/s^2) and mean semimajor axis (km).
    // Only moons with a measured GM are included.
    System mars_system(Body(2, mars_mass, mars_style), o.g);
    add_moon(mars_system, 0.0007087, 9375.0, 0.5, Config{160, 150, 140, 255, 1.108e4f});
    add_moon(mars_system, 0.0000962, 23457.0, 2.4, Config{150, 145, 140, 255, 6.2e3f});

    System jupiter_system(Body(2, jupiter_mass, jupiter_style), o.g);
    add_moon(jupiter_system, 5959.91547, 421800.0, 0.4, Config{230, 210, 80, 255, 1.82149e6f});
    add_moon(jupiter_system, 3202.71210, 671100.0, 1.6, Config{210, 180, 140, 255, 1.56080e6f});
    add_moon(jupiter_system, 9887.83275, 1070400.0, 2.8, Config{170, 160, 140, 255, 2.63120e6f});
    add_moon(jupiter_system, 7179.28340, 1882700.0, 4.2, Config{120, 110, 90, 255, 2.41030e6f});

    System saturn_system(Body(2, saturn_mass, saturn_style), o.g);
    add_moon(saturn_system, 2.50349, 186000.0, 0.3, Config{200, 200, 200, 255, 1.9820e5f});
    add_moon(saturn_system, 7.21037, 238400.0, 1.1, Config{220, 230, 235, 255, 2.5210e5f});
    add_moon(saturn_system, 41.21353, 295000.0, 1.9, Config{210, 210, 205, 255, 5.3110e5f});
    add_moon(saturn_system, 73.11607, 377700.0, 2.7, Config{205, 205, 200, 255, 5.6140e5f});
    add_moon(saturn_system, 153.94175, 527200.0, 3.5, Config{195, 195, 190, 255, 7.6350e5f});
    add_moon(saturn_system, 8978.13710, 1221900.0, 4.4, Config{214, 164, 72, 255, 2.57476e6f});
    add_moon(saturn_system, 120.51511, 3561700.0, 5.4, Config{180, 170, 150, 255, 7.3430e5f});

    System uranus_system(Body(2, uranus_mass, uranus_style), o.g);
    add_moon(uranus_system, 4.3, 129846.0, 0.6, Config{180, 190, 195, 255, 2.358e5f});
    add_moon(uranus_system, 83.5, 190929.0, 1.8, Config{200, 205, 210, 255, 5.789e5f});
    add_moon(uranus_system, 85.1, 265986.0, 3.0, Config{140, 145, 150, 255, 5.847e5f});
    add_moon(uranus_system, 226.9, 436298.0, 4.2, Config{190, 195, 200, 255, 7.889e5f});
    add_moon(uranus_system, 205.3, 583511.0, 5.4, Config{170, 165, 160, 255, 7.614e5f});

    System neptune_system(Body(2, neptune_mass, neptune_style), o.g);
    // Triton is retrograde. The last argument flips the orbit direction.
    add_moon(neptune_system, 1428.49546, 354800.0, 1.2, Config{190, 180, 180, 255, 1.3526e6f}, true);

    // add() places the Earth-Moon center of mass. Shift that center outward
    // by the barycenter offset so the Earth itself, not the barycenter, is
    // one astronomical unit from the Sun. The Moon lies further out on the
    // same line.
    const double barycenter_offset =
        moon_mass / (earth_mass + moon_mass) * earth_moon_distance;

    // Position angles are initial phases, so the planets do not start on
    // one line. Distances are the circular-orbit radii. Earth stays at
    // angle 0 so the Moon remains directly outside it. Each planet with
    // moons is placed by its barycenter.
    System solar_system(Body(2, sun_mass, sun_style), o.g);
    solar_system.add(Body(2, mercury_mass, mercury_style), mercury_distance, 0.7);
    solar_system.add(Body(2, venus_mass, venus_style), venus_distance, 2.0);
    solar_system.add(
        std::move(earth_system),
        earth_sun_distance + barycenter_offset,
        0.0
    );
    solar_system.add(std::move(mars_system), mars_distance, 2.8);
    solar_system.add(std::move(jupiter_system), jupiter_distance, 3.7);
    solar_system.add(std::move(saturn_system), saturn_distance, 4.5);
    solar_system.add(std::move(uranus_system), uranus_distance, 5.3);
    solar_system.add(std::move(neptune_system), neptune_distance, 6.0);

    const State init_state = solar_system.realize();

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
        "Solar System",
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
